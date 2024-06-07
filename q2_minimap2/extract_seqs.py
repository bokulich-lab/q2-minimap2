# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile

from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2._filtering_utils import (
    convert_to_fasta,
    make_mn2_cmd,
    process_sam_file,
    run_cmd,
    set_penalties,
)
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


# This function uses Minimap2 to align reads to a reference, filters the
# resulting SAM file based on mapping criteria using samtools view, and
# finally converts the filtered SAM file to FASTA format using samtools fasta,
# saving the output in the specified directory.
def _minimap2_extract_seqs(
    reads,
    outdir,
    idx_path,
    n_threads,
    preset,
    keep,
    min_per_identity,
    penalties,
):

    with tempfile.NamedTemporaryFile() as sam_f:
        # Use Minimap2 to find mapped and unmapped reads
        mn2_cmd = make_mn2_cmd(
            preset, idx_path, n_threads, penalties, reads, None, sam_f.name
        )
        run_cmd(mn2_cmd, "Minimap2")

        # Filter sam file using samtools view
        process_sam_file(sam_f.name, keep, min_per_identity)

        # Convert to FASTA with samtools
        convert_to_fasta_cmd = convert_to_fasta(outdir.path, n_threads, sam_f.name)
        run_cmd(convert_to_fasta_cmd, "samtools fasta")


def extract_seqs(
    sequences: DNAFASTAFormat,
    index_database: Minimap2IndexDBDirFmt = None,
    reference_reads: DNAFASTAFormat = None,
    n_threads: int = 3,
    mapping_preset: str = "map-ont",
    extract: str = "mapped",
    min_per_identity: float = None,
    matching_score: int = None,
    mismatching_penalty: int = None,
    gap_open_penalty: int = None,
    gap_extension_penalty: int = None,
) -> DNAFASTAFormat:

    if reference_reads and index_database:
        raise ValueError(
            "Only one reference_reads or index_database artifact "
            "can be provided as input. Choose one and try again."
        )

    # Initialize directory format for filtered sequences
    filtered_seqs = DNAFASTAFormat()

    # Get the path of index database or reference
    if index_database:
        idx_ref_path = str(index_database.path) + "/index.mmi"
    elif reference_reads:
        idx_ref_path = str(reference_reads.path)
    else:
        raise ValueError(
            "Either reference_reads or a minimap2_index must be provided as input."
        )

    penalties = set_penalties(
        matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
    )

    # Filter the read using minimap2 according to the specified parameters
    _minimap2_extract_seqs(
        str(sequences),
        filtered_seqs,
        idx_ref_path,
        n_threads,
        mapping_preset,
        extract,
        min_per_identity,
        penalties,
    )

    return filtered_seqs
