# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2._filtering_utils import run_cmd
from q2_minimap2.types._format import Minimap2IndexDBDirFmt, PairwiseAlignmentMN2Format


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
    # Filter mapped query entries based on identity score
    mapped_df = df[df[9] / df[10] >= perc_identity]

    if output_no_hits:
        filtered_out = df[(df[9] / df[10] < perc_identity) | (df[10] == 0)]

        # Keep only the first entry/row for each unique query
        # that is filtered out, which signifies that we treat it
        # as an unmapped query
        filtered_out = filtered_out.drop_duplicates(subset=1)

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
    idx_ref_path, query_reads, n_threads, mapping_preset, paf_file_fp, output_no_hits
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
    n_threads: int = 3,
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

    # Determine the reference or index path based on input
    idx_ref_path = str(index.path / "index.mmi") if index else str(reference.path)

    # Create a reference to a file with PAF format
    paf_file_fp = PairwiseAlignmentMN2Format()

    # Construct the command
    cmd = construct_command(
        idx_ref_path,
        query,
        n_threads,
        preset,
        paf_file_fp,
        output_no_hits,
    )

    # Execute the Minimap2 alignment command
    run_cmd(cmd, "Minimap2")

    # Read the PAF file as a pandas DataFrame
    df = pd.read_csv(str(paf_file_fp), sep="\t", header=None)

    # Filter the PAF file by maxaccepts (default = 1)
    df = filter_by_maxaccepts(df, maxaccepts)

    # Optionally filter by perc_identity
    if min_per_identity is not None:
        df = filter_by_perc_identity(df, min_per_identity, output_no_hits)

    df.reset_index(drop=True, inplace=True)

    return df
