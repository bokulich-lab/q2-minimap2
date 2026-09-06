# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat
from qiime2.plugin import get_available_cores

from q2_minimap2._filtering_utils import run_cmd
from q2_minimap2.types._format import (
    Minimap2IndexDBDirFmt,
    PairwiseAlignmentMN2Format,
    read_paf,
)


# Filter a PAF file to keep only a certain number of entries
# for each individual query name, as defined by "maxaccepts"
def filter_by_maxaccepts(df, maxaccepts):
    # Group by query_name and count occurrences
    counts = df.groupby(0).cumcount() + 1

    # Filter the DataFrame based on maxaccepts
    filtered_df = df[counts <= maxaccepts]

    return filtered_df


# Filter PAF entries based on a threshold of percentage identity
def filter_by_perc_identity(df, perc_identity, output_no_hits):
    # Filter mapped query entries based on identity score. No-hit rows have a
    # block length of 0, so this division yields NaN for them and the
    # comparison is False, which keeps them out of the mapped set.
    mapped_df = df[df[9] / df[10] >= perc_identity]

    if output_no_hits:
        # Only queries left without a single accepted alignment become no-hit
        # rows. Selecting on the query name (column 0) keeps exactly one
        # placeholder per query; de-duplicating on any other column would drop
        # distinct queries that happen to share a value.
        filtered_out = df[~df[0].isin(mapped_df[0])]
        filtered_out = filtered_out.drop_duplicates(subset=0).copy()

        # The strand and target name columns hold "*" in a no-hit row, so they
        # have to accept strings even when every reference ID is numeric and
        # pandas typed the column as an integer.
        filtered_out[[4, 5]] = filtered_out[[4, 5]].astype(object)

        # Change paf file column entries that are filtered out
        # to indicate that are unmapped queries
        filtered_out.iloc[:, 2:12] = 0
        filtered_out.iloc[:, 4:6] = "*"
        filtered_out.iloc[:, 12:] = "*"

        # Merging the two DataFrames based on row number
        mapped_df = pd.concat([mapped_df, filtered_out], axis=0)
        mapped_df = mapped_df.sort_index()

    return mapped_df


# Construct the command list for the Minimap2 alignment search
def construct_command(
    idx_ref_path,
    query_reads,
    n_threads,
    mapping_preset,
    paf_file_fp,
    output_no_hits,
    maxaccepts=1,
    split_prefix=None,
):
    cmd = [
        "minimap2",
        "-x",
        mapping_preset,
        "-c",
        str(idx_ref_path),
        str(query_reads),
        "-t",
        str(n_threads),
        # Minimap2 retains only 5 secondary alignments by default, which caps
        # every query at 6 hits and makes any larger maxaccepts unreachable.
        # Asking for maxaccepts secondaries leaves enough rows to filter down to.
        "-N",
        str(maxaccepts),
        # Needed so a multi-part index still yields globally ranked hits rather
        # than hits grouped per index part. No-op for a single-part index.
        "--split-prefix",
        str(split_prefix),
        "-o",
        str(paf_file_fp),
    ]
    if output_no_hits:
        cmd.append("--paf-no-hit")
    return cmd


# Performs sequence alignment using Minimap2 and outputs results in
# PairwiseAlignmentMN2 format.
def minimap2_search(
    query: DNAFASTAFormat,
    index: Minimap2IndexDBDirFmt = None,
    reference: DNAFASTAFormat = None,
    n_threads: int = 1,
    preset: str = "map-ont",
    maxaccepts: int = 1,
    min_per_identity: float = None,
    output_no_hits: bool = True,
) -> pd.DataFrame:
    # Ensure that only one of reference or index is provided
    if reference and index:
        raise ValueError(
            "Only one of reference or index can be provided as input. "
            "Choose one and try again."
        )

    # Ensure that at least one of reference and index is provided
    if not reference and not index:
        raise ValueError("Either reference or index must be provided as input.")

    # A thread count of 0 is what the Threads type passes on for "auto"
    if n_threads == 0:
        n_threads = get_available_cores()

    # Determine the reference or index path based on input
    idx_ref_path = str(index.path / "index.mmi") if index else str(reference.path)

    # Create a reference to a file with PAF format
    paf_file_fp = PairwiseAlignmentMN2Format()

    with tempfile.TemporaryDirectory() as tmpd:
        # Construct the command
        cmd = construct_command(
            idx_ref_path,
            query,
            n_threads,
            preset,
            paf_file_fp,
            output_no_hits,
            maxaccepts,
            os.path.join(tmpd, "split"),
        )

        # Execute the Minimap2 alignment command
        run_cmd(cmd, "Minimap2")

    # Read the PAF file as a pandas DataFrame
    df = read_paf(str(paf_file_fp))

    # Optionally filter by perc_identity. This runs before the maxaccepts cut
    # so that maxaccepts selects among the hits that actually pass the identity
    # threshold, which is what the parameter description promises. Truncating
    # first would discard qualifying hits in favour of higher-scoring ones that
    # the threshold then rejects.
    if min_per_identity is not None:
        df = filter_by_perc_identity(df, min_per_identity, output_no_hits)

    # Filter the PAF file by maxaccepts (default = 1)
    df = filter_by_maxaccepts(df, maxaccepts)

    df.reset_index(drop=True, inplace=True)

    return df
