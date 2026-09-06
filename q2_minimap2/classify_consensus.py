# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from q2_feature_classifier._consensus_assignment import _compute_consensus_annotations


def _PAF_format_df_to_series_of_lists(
    assignments: pd.DataFrame,  # Input dataframe containing PAF formatted data
    ref_taxa: pd.Series,  # Series mapping sequence IDs to taxonomic labels
    unassignable_label: str = "Unassigned",  # Label for indeterminate assignments
) -> pd.Series:

    # Only the query and subject identifier columns take part in the lookup, so
    # converting those is enough; converting the whole frame would copy the
    # CIGAR column, which is by far the largest, for nothing
    query_ids = assignments[0].astype(str)
    subject_ids = assignments[5].astype(str)

    # Identify missing IDs in the reference taxonomy that are present in the
    # assignments, excluding special cases ('*' for unassigned and empty strings).
    missing_ids = set(subject_ids.values) - set(ref_taxa.index) - {"*", ""}
    if len(missing_ids) > 0:
        # Raise an error if there are IDs found in assignments that aren't in the
        # reference taxonomy
        raise KeyError(
            "Reference taxonomy and search results do not match. "
            "The following identifiers were reported in the search "
            "results but are not present in the reference taxonomy:"
            " {0}".format(", ".join(str(i) for i in missing_ids))
        )

    # Assign a label to '*' (representing sequences without assignment) in the
    # reference taxonomy. This works on a copy so that the caller's taxonomy is
    # left as it was.
    ref_taxa = ref_taxa.copy()
    ref_taxa["*"] = unassignable_label

    # Replace every subject ID with its taxonomic label in one pass. Doing this
    # row by row also meant indexing positionally with a label, which silently
    # mis-assigned taxonomy for any caller whose frame was not indexed from 0.
    taxa_hits: pd.Series = subject_ids.map(ref_taxa)
    taxa_hits.index = query_ids

    # Group all annotations for each unique accession ID together
    taxa_hits = taxa_hits.groupby(taxa_hits.index).apply(list)

    return taxa_hits


def _find_consensus_annotation(
    search_results: pd.DataFrame,  # PAF search results
    reference_taxonomy: pd.Series,  # Mapping of sequence IDs to taxonomy
    min_consensus: int = 0.51,  # Threshold for consensus decision
    unassignable_label: str = "Unassigned",  # Label for indeterminate assignments
) -> pd.DataFrame:

    # Convert the PAF results to a series mapping each feature to a list of observed
    # taxonomic hits. This is a preprocessing step that formats the search results
    # for consensus computation
    obs_taxa = _PAF_format_df_to_series_of_lists(
        search_results, reference_taxonomy, unassignable_label=unassignable_label
    )

    # Compute the consensus annotations for each feature based on the observed
    # taxonomic hits. This involves finding the most common taxonomy for each
    # feature and comparing it against the minimum consensus threshold to decide
    # if a consensus taxonomy can be assigned
    result = _compute_consensus_annotations(
        obs_taxa, min_consensus=min_consensus, unassignable_label=unassignable_label
    )

    # Set the name of the DataFrame index to "Feature ID", indicating that each row
    # corresponds to a unique feature for which a consensus taxonomy has been
    # determined (or marked unassignable)
    result.index.name = "Feature ID"

    return result


def classify_consensus_minimap2(
    ctx,
    query,  # Query sequences for classification
    reference_taxonomy,  # Taxonomy mappings for reference sequences
    index=None,  # Optional pre-built index of the reference database
    reference=None,  # Optional reference sequences for on-the-fly indexing
    maxaccepts=10,  # Maximum number of alignments to accept per query
    preset="map-ont",  # Minimap2 mapping preset
    min_per_identity=0.8,  # Min percentage identity for an alignment to be accepted
    output_no_hits=True,  # Include queries with no hits in the output
    min_consensus=0.51,  # Threshold for consensus in classification
    unassignable_label="Unassigned",  # Label for unclassifiable queries
    n_threads=3,  # Number of threads for the alignment process
):
    # Retrieve the necessary actions from the context to perform
    # sequence mapping and consensus annotation
    search_db = ctx.get_action("minimap2", "minimap2_search")
    lca = ctx.get_action("minimap2", "_find_consensus_annotation")

    # Execute the search against a reference database or reads using Minimap2,
    # adhering to specified parameters like maxaccepts and perc_identity
    (result,) = search_db(
        query=query,
        index=index,
        reference=reference,
        maxaccepts=maxaccepts,
        preset=preset,
        min_per_identity=min_per_identity,
        output_no_hits=output_no_hits,
        n_threads=n_threads,
    )

    # Utilize the search results to find consensus annotations for each query sequence,
    # leveraging the reference taxonomy and consensus threshold provided
    (consensus,) = lca(
        search_results=result,
        reference_taxonomy=reference_taxonomy,
        min_consensus=min_consensus,
        unassignable_label=unassignable_label,
    )

    return result, consensus
