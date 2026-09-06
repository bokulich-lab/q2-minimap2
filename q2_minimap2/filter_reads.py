# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2._filtering_utils import (
    collate_sam_inplace,
    convert_to_fastq,
    make_mn2_cmd,
    process_paired_sam_flags,
    process_sam_file,
    run_cmd,
    set_penalties,
)
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


# This function aligns reads to a reference using Minimap2,
# filters the alignments, and converts the filtered alignments
# to FASTQ format, storing the output in a given directory
def _minimap2_filter_reads(
    reads1,
    reads2,
    outdir,  # Output directory for filtered reads in FASTQ format
    idx_path,  # Path to the Minimap2 index file
    n_threads,  # Number of threads for Minimap2 and samtools
    preset,  # Preset options for Minimap2 alignment
    keep,  # Flag to exclude mapped reads, keeping only unmapped
    min_per_identity,  # Minimum percentage identity for a read to be included
    penalties,  # Scoring penalties for Minimap2 alignment
):
    # Create a temporary file for SAM output from Minimap2, plus a scratch
    # directory for the prefix Minimap2 needs to merge a multi-part index
    with tempfile.NamedTemporaryFile() as sam_f, tempfile.TemporaryDirectory() as tmpd:
        # Construct and execute Minimap2 command for alignment
        mn2_cmd = make_mn2_cmd(
            preset,
            idx_path,
            n_threads,
            penalties,
            reads1,
            reads2,
            sam_f.name,
            os.path.join(tmpd, "split"),
        )

        run_cmd(mn2_cmd, "Minimap2")

        # Filter the SAM file using samtools based on the include/exclude criteria
        process_sam_file(sam_f.name, keep, min_per_identity)

        if reads2:
            # Ensuring proper read grouping of paired reads
            collate_sam_inplace(sam_f.name)
            # Making flags suitable for samtools fastq command
            process_paired_sam_flags(sam_f.name)

            # Construct and execute command to convert SAM to FASTQ
            # using samtools fastq, directing output to the specified output directory
            file1 = str(outdir.path / os.path.basename(reads1))
            file2 = str(outdir.path / os.path.basename(reads2))
            _reads = ["-1", file1, "-2", file2]
            convert_to_fastq_cmd = convert_to_fastq(
                _reads, n_threads, sam_f.name, "paired"
            )
        else:
            # Convert SAM to FASTQ
            fwd = str(outdir.path / os.path.basename(reads1))
            _reads = ["-0", fwd]
            convert_to_fastq_cmd = convert_to_fastq(
                _reads, n_threads, sam_f.name, "single"
            )

        run_cmd(convert_to_fastq_cmd, "samtools fastq")


def filter_reads(
    query: CasavaOneEightSingleLanePerSampleDirFmt,
    index: Minimap2IndexDBDirFmt = None,  # Optional pre-built Minimap2 index
    reference: DNAFASTAFormat = None,  # Optional reference sequences
    n_threads: int = 3,  # Number of threads for Minimap2
    preset: str = "map-ont",  # Minimap2 mapping preset
    keep: str = "mapped",  # Keep 'mapped' or 'unmapped' reads
    min_per_identity: float = None,  # Minimum percentage identity to keep a read
    matching_score: int = None,  # Score for matching bases
    mismatching_penalty: int = None,  # Penalty for mismatched bases
    gap_open_penalty: int = None,  # Penalty for opening a gap
    gap_extension_penalty: int = None,  # Penalty for extending a gap
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    # Ensure that only one of reference or index is provided
    if reference and index:
        raise ValueError(
            "Only one of reference or index can be provided as input. "
            "Choose one and try again."
        )

    # Ensure that at least one of reference and index is provided
    if not reference and not index:
        raise ValueError("Either reference or index must be provided as input.")

    # Determine the reference path from the provided index or reference
    if index:
        idx_ref_path = str(index) + "/index.mmi"  # Path to index file
    elif reference:
        idx_ref_path = str(reference)  # Path to the reference file

    # Initialize directory format for filtered sequences
    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()

    # Set penalties for alignment based on function parameters
    penalties = set_penalties(
        matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
    )

    # Process each read, filtering according to the specified parameters
    for _, fwd, rev in query.manifest.itertuples():
        _minimap2_filter_reads(
            fwd,
            rev,
            filtered_seqs,
            idx_ref_path,
            n_threads,
            preset,
            keep,
            min_per_identity,
            penalties,
        )

    return filtered_seqs
