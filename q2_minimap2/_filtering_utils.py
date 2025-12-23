# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import re
import shutil
import subprocess
import tempfile

from q2_minimap2._utils import run_command


# Set Minimap2 alignment penalties based on provided parameters
def set_penalties(
    matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
):
    options = []
    if matching_score is not None:
        options += ["-A", str(matching_score)]
    if mismatching_penalty is not None:
        options += ["-B", str(mismatching_penalty)]
    if gap_open_penalty:
        options += ["-O", str(gap_open_penalty)]
    if gap_extension_penalty:
        options += ["-E", str(gap_extension_penalty)]

    return options


# Function to calculate the identity percentage of an alignment
def calculate_identity(aln, total_length):
    try:
        # Extracts the number of mismatches (NM tag) from the SAM file alignment line
        nm = int([x for x in aln.split("\t") if x.startswith("NM:i:")][0].split(":")[2])
    except IndexError:
        # Defaults to 0 mismatches if the NM tag is not found
        nm = 0

    # Calculates matches by subtracting mismatches from total length
    matches = total_length - nm

    # Calculates identity percentage as the ratio of matches to total alignment length.
    identity_percentage = matches / total_length

    return identity_percentage


# Function to get the alignment length from a CIGAR string
def get_alignment_length(cigar):
    if cigar == "*":
        # Returns 0 if the CIGAR string is '*', indicating no alignment
        return 0

    # Extracts all match, insertion, and deletion operations from the CIGAR string
    matches = re.findall(r"(\d+)([MID])", cigar)

    # Sums the lengths of matches, insertions, and deletions to get total
    # alignment length
    total_length = sum(int(length) for length, op in matches if op in ["M", "D", "I"])

    return total_length


# Function to process a SAM file, filter based on mappings and identity percentage
def process_sam_file(input_sam_file, keep, min_per_identity):
    # Creates a temporary file and opens the input SAM file for reading simultaneously
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file, open(
        input_sam_file, "r"
    ) as infile:
        for line in infile:
            # Writes header lines directly to the output file
            if line.startswith("@"):
                tmp_file.write(line)
                continue

            # Extract information from the line
            parts = line.split("\t")
            flag = int(parts[1])
            cigar = parts[5]

            # Calculates identity percentage for alignments with a valid CIGAR string
            if min_per_identity and cigar != "*":
                total_length = get_alignment_length(cigar)
                identity_percentage = calculate_identity(line, total_length)
            else:
                # Defaults identity percentage to 100% if no CIGAR string or no
                # min_per_identity specified
                identity_percentage = 1

            # Logic for including or excluding reads based on mappings and
            # identity percentage
            if keep == "mapped":
                if not (flag & 0x4) and not (flag & 0x100):
                    if not min_per_identity or identity_percentage >= min_per_identity:
                        tmp_file.write(line)
            else:
                # Condition for keeping unmapped reads or mapped reads below the
                # identity threshold
                if (flag & 0x4) or (
                    min_per_identity and identity_percentage < min_per_identity
                ):
                    tmp_file.write(line)

    # Replaces the original SAM file with the filtered temporary file
    shutil.move(tmp_file.name, input_sam_file)


# Generate samtools fasta convert command
def convert_to_fasta(_reads, n_threads, samfile_filepath):
    # -s /dev/null excludes singletons
    # -n keeps samtools from altering header IDs
    convert_cmd = [
        "samtools",
        "fasta",
        "-0",
        str(_reads),
        "-s",
        "/dev/null",
        "-@",
        str(n_threads),
        "-n",
        str(samfile_filepath),
    ]

    return convert_cmd


def convert_to_fastq(_reads, n_threads, samfile_filepath, kind):
    convert_cmd = ["samtools", "fastq", *_reads]
    if kind == "paired":
        convert_cmd += ["-0", "/dev/null"]

    convert_cmd += [
        "-s",
        "/dev/null",
        "-@",
        str(n_threads),
        "-n",
        str(samfile_filepath),
    ]

    return convert_cmd


# Generate Minimap2 mapping command
def make_mn2_cmd(mapping_preset, index, n_threads, penalties, reads1, reads2, samf_fp):
    # align to reference with Minimap2
    minimap2_cmd = (
        [
            "minimap2",
            "-a",
            "-x",
            mapping_preset,
            str(index),
            "-t",
            str(n_threads),
            "-o",
            str(samf_fp),
        ]
        + penalties
        + [reads1]
    )

    if reads2:
        minimap2_cmd.append(reads2)

    return minimap2_cmd


# Helper function for command execution
def run_cmd(cmd, tool_name):
    try:
        # Execute samtools fastq
        run_command(cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while using {tool_name}, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def collate_sam_inplace(input_sam_path):
    # Temporary file prefix based on the input file name
    temp_prefix = os.path.splitext(input_sam_path)[0] + "_temp_collate"
    # Output file path in the same directory
    output_sam_path = os.path.splitext(input_sam_path)[0] + "_collated.sam"

    # Samtools collate command
    collate_cmd = [
        "samtools",
        "collate",
        "-u",
        "-o",
        output_sam_path,
        "-T",
        temp_prefix,
        input_sam_path,
    ]

    # Execute samtools collate command
    run_cmd(collate_cmd, "Samtools collate")

    shutil.move(output_sam_path, input_sam_path)


def process_paired_sam_flags(input_sam_path):
    """
    Process a SAM file containing paired-end reads to set specific flags for the read
    pairs in order to be recognized py the samtools fastq command for paired end reads
    """
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as temp_file:
        with open(input_sam_path, "r") as infile:
            for line in infile:
                if line.startswith("@"):
                    temp_file.write(line)
                    continue

                read1 = line.strip().split("\t")
                read2 = infile.readline().strip().split("\t")

                if int(read1[1]) == 4 and int(read2[1]) == 4:  # Both reads are unmapped
                    read1[1] = "69"  # 1 + 4 + 64: Read is first in pair and unmapped
                    read2[1] = "133"  # 1 + 4 + 128: Read is second in pair and unmapped
                else:
                    # Ensure paired flags
                    # Add 1 and 64 (first read in a pair)
                    read1[1] = str(int(read1[1]) | 1 | 64)
                    # Add 1 and 128 (first read in a pair)
                    read2[1] = str(int(read2[1]) | 1 | 128)

                temp_file.write("\t".join(read1) + "\n")
                temp_file.write("\t".join(read2) + "\n")

    shutil.move(temp_file.name, input_sam_path)
