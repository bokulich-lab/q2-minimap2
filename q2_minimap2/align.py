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
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    CasavaOneEightSingleLanePerSampleDirFmt,
)
from qiime2.plugin import get_available_cores

from q2_minimap2._filtering_utils import (
    get_index_or_reference_path,
    make_mn2_cmd,
    run_cmd,
    set_penalties,
)
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


# Align one sample and store the result as a coordinate-sorted BAM named after
# the sample. The SAM is only ever an intermediate: BAM keeps the alignment
# compressed, which matters because a long-read SAM is very large.
def _minimap2_align_sample(
    reads1,
    reads2,
    sample_id,
    outdir,
    idx_path,
    n_threads,
    preset,
    penalties,
):
    with tempfile.NamedTemporaryFile() as sam_f, tempfile.TemporaryDirectory() as tmpd:
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

        sort_cmd = [
            "samtools",
            "sort",
            "-@",
            str(n_threads),
            "-o",
            str(outdir.path / f"{sample_id}.bam"),
            sam_f.name,
        ]
        run_cmd(sort_cmd, "samtools sort")


def align(
    query: CasavaOneEightSingleLanePerSampleDirFmt,
    index: Minimap2IndexDBDirFmt = None,
    reference: DNAFASTAFormat = None,
    n_threads: int = 1,
    preset: str = "map-ont",
    matching_score: int = None,
    mismatching_penalty: int = None,
    gap_open_penalty: int = None,
    gap_extension_penalty: int = None,
) -> BAMDirFmt:

    idx_ref_path = get_index_or_reference_path(index, reference)

    # A thread count of 0 is what the Threads type passes on for "auto"
    if n_threads == 0:
        n_threads = get_available_cores()

    alignments = BAMDirFmt()

    penalties = set_penalties(
        matching_score, mismatching_penalty, gap_open_penalty, gap_extension_penalty
    )

    for sample_id, fwd, rev in query.manifest.itertuples():
        _minimap2_align_sample(
            fwd,
            rev,
            sample_id,
            alignments,
            idx_ref_path,
            n_threads,
            preset,
            penalties,
        )

    return alignments
