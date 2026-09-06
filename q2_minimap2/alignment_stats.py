# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

import pandas as pd
import qiime2
from q2_types.per_sample_sequences import BAMDirFmt


# Run a samtools subcommand and hand back its stdout. The alignment summaries
# are small, so there is nothing to gain from streaming them through a file.
def _run_samtools(args):
    result = subprocess.run(args, check=True, capture_output=True, text=True)

    return result.stdout


# Parse "samtools flagstat -O tsv", whose rows are
# <QC-passed count>\t<QC-failed count>\t<metric name>
def _parse_flagstat(output):
    counts = {}
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue

        try:
            counts[fields[2]] = int(fields[0])
        except ValueError:
            # The percentage rows are reported as "83.33%" or "N/A", and both
            # are derivable from the counts, so they are of no use here
            continue

    return counts


# Parse "samtools coverage", one row per reference sequence
def _parse_coverage(output):
    rows = []
    for line in output.splitlines():
        if line.startswith("#") or not line.strip():
            continue

        fields = line.split("\t")
        rows.append(
            {
                "numreads": int(fields[3]),
                "coverage": float(fields[5]),
                "meandepth": float(fields[6]),
                "meanmapq": float(fields[8]),
            }
        )

    return rows


# Collect the per-sample summary for one BAM
def _summarize_bam(bam_fp):
    counts = _parse_flagstat(
        _run_samtools(["samtools", "flagstat", "-O", "tsv", bam_fp])
    )
    references = _parse_coverage(_run_samtools(["samtools", "coverage", bam_fp]))

    total = counts.get("primary", 0)
    mapped = counts.get("primary mapped", 0)
    with_reads = [ref for ref in references if ref["numreads"] > 0]

    # Mapping quality is averaged over the reads rather than over the
    # references, so that a reference carrying a single read does not weigh as
    # much as one carrying thousands
    mapped_in_refs = sum(ref["numreads"] for ref in with_reads)
    if mapped_in_refs > 0:
        mean_mapq = (
            sum(ref["meanmapq"] * ref["numreads"] for ref in with_reads)
            / mapped_in_refs
        )
    else:
        mean_mapq = 0.0

    return {
        "total_reads": total,
        "mapped_reads": mapped,
        "unmapped_reads": total - mapped,
        "percent_mapped": (100 * mapped / total) if total else 0.0,
        "secondary_alignments": counts.get("secondary", 0),
        "supplementary_alignments": counts.get("supplementary", 0),
        "references_total": len(references),
        "references_covered": len(with_reads),
        "mean_mapping_quality": mean_mapq,
        "mean_depth": (
            sum(ref["meandepth"] for ref in references) / len(references)
            if references
            else 0.0
        ),
        "percent_reference_covered": (
            sum(ref["coverage"] for ref in references) / len(references)
            if references
            else 0.0
        ),
    }


def alignment_stats(alignments: BAMDirFmt) -> qiime2.Metadata:
    summaries = {}
    for bam_fp in sorted(alignments.path.glob("*.bam")):
        summaries[bam_fp.stem] = _summarize_bam(str(bam_fp))

    stats = pd.DataFrame.from_dict(summaries, orient="index").astype(float)
    stats.index.name = "sample-id"

    return qiime2.Metadata(stats)
