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

# SAM FLAG bits used when filtering and re-pairing alignment records
PAIRED = 0x1
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
FIRST_IN_PAIR = 0x40
SECOND_IN_PAIR = 0x80
SECONDARY = 0x100
SUPPLEMENTARY = 0x800


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
    tmp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    try:
        with tmp_file, open(input_sam_file, "r") as infile:
            for line in infile:
                # Writes header lines directly to the output file
                if line.startswith("@"):
                    tmp_file.write(line)
                    continue

                # Extract information from the line
                parts = line.split("\t")
                flag = int(parts[1])
                cigar = parts[5]

                # Identity percentage, for a valid CIGAR string
                if min_per_identity and cigar != "*":
                    total_length = get_alignment_length(cigar)
                    identity_percentage = calculate_identity(line, total_length)
                else:
                    # Defaults identity percentage to 100% if no CIGAR string or no
                    # min_per_identity specified
                    identity_percentage = 1

                # A secondary or supplementary record is an extra alignment of a
                # read that its primary record already represents, so keeping it
                # would emit the same read more than once
                if flag & SECONDARY or flag & SUPPLEMENTARY:
                    continue

                # Logic for including or excluding reads based on mappings and
                # identity percentage
                if keep == "mapped":
                    if not (flag & UNMAPPED):
                        if (
                            not min_per_identity
                            or identity_percentage >= min_per_identity
                        ):
                            tmp_file.write(line)
                else:
                    # Condition for keeping unmapped reads or mapped reads below the
                    # identity threshold
                    if (flag & UNMAPPED) or (
                        min_per_identity and identity_percentage < min_per_identity
                    ):
                        tmp_file.write(line)

        # Replaces the original SAM file with the filtered temporary file
        shutil.move(tmp_file.name, input_sam_file)
    finally:
        # A run that fails part-way through would otherwise leave behind a
        # temporary copy of the SAM, which for a long-read run is as large as
        # the input itself
        if os.path.exists(tmp_file.name):
            os.remove(tmp_file.name)


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
def make_mn2_cmd(
    mapping_preset, index, n_threads, penalties, reads1, reads2, samf_fp, split_prefix
):
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
    else:
        # A reference larger than the -I threshold produces a multi-part index,
        # and without this Minimap2 emits no @SQ header lines and repeats every
        # read once per part, which makes the samtools steps fail. Verified to
        # leave single-part output byte-identical, so it is safe to pass always
        # here. It is deliberately not passed for paired input: there it also
        # switches Minimap2 into fragment mode, which changes which reads align.
        minimap2_cmd += ["--split-prefix", str(split_prefix)]

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


# Rewrite the flags of a single mate pair so that the samtools fastq command
# recognises it. Minimap2 usually sets the first/second-in-pair bits itself, in
# which case they are kept; they are only assigned from file order when absent.
def set_pair_flags(read1, read2):
    flag1, flag2 = int(read1[1]), int(read2[1])

    # Clear the mate bits before setting them, so that a record which already
    # carries one of them cannot end up flagged as both first and second in
    # pair, which samtools fastq discards.
    base1 = flag1 & ~(PAIRED | FIRST_IN_PAIR | SECOND_IN_PAIR | MATE_UNMAPPED)
    base2 = flag2 & ~(PAIRED | FIRST_IN_PAIR | SECOND_IN_PAIR | MATE_UNMAPPED)

    new1 = base1 | PAIRED | FIRST_IN_PAIR
    new2 = base2 | PAIRED | SECOND_IN_PAIR

    # Each mate records whether the other one is unmapped
    if flag2 & UNMAPPED:
        new1 |= MATE_UNMAPPED
    if flag1 & UNMAPPED:
        new2 |= MATE_UNMAPPED

    read1[1], read2[1] = str(new1), str(new2)

    return read1, read2


# Order the two records of a pair, preferring the first/second-in-pair bits
# Minimap2 set and falling back to the order they appear in the file.
def order_mates(group):
    first = [read for read in group if int(read[1]) & FIRST_IN_PAIR]
    second = [read for read in group if int(read[1]) & SECOND_IN_PAIR]

    # A record carrying both bits satisfies each selection, which would return
    # it twice and silently drop its partner, so the two must be distinct
    if len(first) == 1 and len(second) == 1 and first[0] is not second[0]:
        return first[0], second[0]

    return group[0], group[1]


# Write out one group of records sharing a read name. Anything that is not
# exactly a pair is dropped: a mate whose partner was removed by filtering has
# nothing left for the samtools fastq command to pair it with.
def write_mate_pair(temp_file, group):
    if len(group) != 2:
        return

    read1, read2 = set_pair_flags(*order_mates(group))
    temp_file.write("\t".join(read1) + "\n")
    temp_file.write("\t".join(read2) + "\n")


def process_paired_sam_flags(input_sam_path):
    """
    Process a SAM file containing paired-end reads to set specific flags for the read
    pairs in order to be recognized py the samtools fastq command for paired end reads.

    The file has already been name-collated, so all the records of a read name are
    consecutive. Grouping by name instead of consuming two lines at a time keeps the
    pairing correct when filtering has removed one mate, or when a read produced more
    than the two records the pair consists of.
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    try:
        with temp_file, open(input_sam_path, "r") as infile:
            group, group_name = [], None

            for line in infile:
                if line.startswith("@"):
                    temp_file.write(line)
                    continue

                read = line.rstrip("\n").split("\t")
                flag = int(read[1])

                # A secondary or supplementary record is an additional alignment
                # of a read that its primary record already stands for, so
                # keeping it here would break the pairing
                if flag & SECONDARY or flag & SUPPLEMENTARY:
                    continue

                if read[0] != group_name:
                    write_mate_pair(temp_file, group)
                    group, group_name = [], read[0]

                group.append(read)

            write_mate_pair(temp_file, group)

        shutil.move(temp_file.name, input_sam_path)
    finally:
        # Do not leave a full copy of the SAM behind if this fails part-way
        if os.path.exists(temp_file.name):
            os.remove(temp_file.name)
