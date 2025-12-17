# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from q2_types.per_sample_sequences import (
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Bool, Choices, Float, Int, Range, Str, TypeMatch

from q2_minimap2.types._type import Minimap2IndexDB, PairwiseAlignmentMN2

T = TypeMatch([SequencesWithQuality, PairedEndSequencesWithQuality])

# filter_reads
filter_reads_inputs = {
    "query": SampleData[T],
    "index": Minimap2IndexDB,
    "reference": FeatureData[Sequence],
}
filter_reads_outputs = [("filtered_query", SampleData[T])]
filter_reads_inputs_dsc = {
    "query": "Sequences to be filtered.",
    "index": "Minimap2 index database. Incompatible with reference.",
    "reference": "Reference sequences. Incompatible with index.",
}
filter_reads_outputs_dsc = {
    "filtered_query": "The resulting filtered sequences.",
}
filter_reads_params = {
    "n_threads": Int % Range(1, None),
    "preset": Str % Choices(["map-ont", "map-hifi", "map-pb", "sr"]),
    "keep": Str % Choices(["mapped", "unmapped"]),
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "matching_score": Int,
    "mismatching_penalty": Int,
    "gap_open_penalty": Int % Range(1, None),
    "gap_extension_penalty": Int % Range(1, None),
}
filter_reads_param_dsc = {
    "n_threads": "Number of threads to use.",
    "preset": "The preset parameter applies multiple options at the same time "
    "during the mapping process of Minimap2. 1) map-ont: Align noisy long reads "
    "of ~10% error rate to a reference genome. 2) map-hifi: Align PacBio "
    "high-fidelity (HiFi) reads to a reference genome. 3) map-pb: Align "
    "older PacBio continuous long reads (CLR) to a reference genome. "
    "4) sr: Align short single-end reads.",
    "keep": "Keep the sequences that align to reference. When "
    "set to unmapped it keeps sequences that do not align to the reference "
    "database.",
    "min_per_identity": "After the alignment step, mapped reads will be "
    "reclassified as unmapped if their identity percentage falls below this "
    "value. If not set, there is no reclassification.",
    "matching_score": "Matching score.",
    "mismatching_penalty": "Mismatching penalty.",
    "gap_open_penalty": "Gap open penalty.",
    "gap_extension_penalty": "Gap extension penalty.",
}
filter_reads_dsc = (
    "This method aligns long-read sequencing data (from a FASTQ file) to a set "
    "of reference sequences, identifying sequences that match or do not match "
    "the reference within a specified identity percentage. The alignment is "
    "performed using Minimap2, and the results are processed using Samtools."
)

# extract-reads
extract_reads_inputs = {
    "sequences": FeatureData[Sequence],
    "index": Minimap2IndexDB,
    "reference": FeatureData[Sequence],
}
extract_reads_inputs_dsc = {
    "sequences": "Sequences to be filtered.",
    "index": "Minimap2 index database. Incompatible with reference.",
    "reference": "Reference sequences. Incompatible with index.",
}
extract_reads_outputs = [("extracted_reads", FeatureData[Sequence])]
extract_reads_outputs_dsc = {
    "extracted_reads": "Subset of sequences that are extracted.",
}
extract_reads_params = {
    "n_threads": Int % Range(1, None),
    "preset": Str % Choices(["map-ont", "map-hifi", "map-pb", "sr"]),
    "extract": Str % Choices(["mapped", "unmapped"]),
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "matching_score": Int,
    "mismatching_penalty": Int,
    "gap_open_penalty": Int % Range(1, None),
    "gap_extension_penalty": Int % Range(1, None),
}
extract_reads_param_dsc = {
    "n_threads": "Number of threads to use.",
    "preset": "The preset parameter applies multiple options at the same time "
    "during the mapping process of Minimap2. 1) map-ont: Align noisy long reads "
    "of ~10% error rate to a reference genome. 2) map-hifi: Align PacBio "
    "high-fidelity (HiFi) reads to a reference genome. 3) map-pb: Align "
    "older PacBio continuous long reads (CLR) to a reference genome. "
    "4) sr: Align short single-end reads.",
    "extract": "Extract sequences that map to reference. When "
    "set to unmapped it extracts sequences that do not map to the reference "
    "database.",
    "min_per_identity": "After the alignment step, mapped reads will be "
    "reclassified as unmapped if their identity percentage falls below this "
    "value. If not set, there is no reclassification.",
    "matching_score": "Matching score.",
    "mismatching_penalty": "Mismatching penalty.",
    "gap_open_penalty": "Gap open penalty.",
    "gap_extension_penalty": "Gap extension penalty.",
}
extract_reads_dsc = (
    "This method aligns long-read sequencing data (from a FASTA file) to a set of "
    "reference sequences, identifying sequences that match or do not match the "
    "reference within a specified identity percentage. The alignment is performed "
    "using Minimap2, and the results are processed using Samtools."
)

# build-index
build_index_inputs = {"reference": FeatureData[Sequence]}
build_index_outputs = [("index", Minimap2IndexDB)]
build_index_inputs_dsc = {"reference": "Reference sequences."}
build_index_outputs_dsc = {"index": "Minimap2 index database."}
build_index_params = {
    "preset": Str % Choices(["map-ont", "map-hifi", "map-pb", "sr"]),
}
build_index_param_dsc = {
    "preset": "This option applies multiple settings at the same time during "
    "the indexing process. This value should match the mapping preset value "
    "that is intended to be used in other actions utilizing the created index. "
    "The available presets are: "
    "1) map-ont: Align noisy long reads of ~10% error rate to a reference genome. "
    "2) map-hifi: Align PacBio high-fidelity (HiFi) reads to a reference genome. "
    "3) map-pb: Align older PacBio continuous long (CLR) reads to a reference genome. "
    "4) sr: Align short single-end reads.",
}
build_index_dsc = "Build a Minimap2 index database from reference sequences."


# minimap2-search
minimap2_search_inputs = {
    "query": FeatureData[Sequence],
    "index": Minimap2IndexDB,
    "reference": FeatureData[Sequence],
}
minimap2_search_outputs = [("search_results", FeatureData[PairwiseAlignmentMN2])]
minimap2_search_inputs_dsc = {
    "query": "Query sequences.",
    "index": "Minimap2 index database. Incompatible with reference.",
    "reference": "Reference sequences. Incompatible with index.",
}
minimap2_search_outputs_dsc = {
    "search_results": "Top hits for each query.",
}
minimap2_search_param_dsc = {
    "n_threads": "Number of threads to use.",
    "maxaccepts": "Maximum number of hits to keep for each query. Minimap2 will "
    "choose the first N hits in the reference database.",
    "preset": "The preset parameter applies multiple options at the same time "
    "during the mapping process of Minimap2. 1) map-ont: Align noisy long reads "
    "of ~10% error rate to a reference genome. 2) map-hifi: Align PacBio "
    "high-fidelity (HiFi) reads to a reference genome. 3) map-pb: Align "
    "older PacBio continuous long reads (CLR) to a reference genome. "
    "4) sr: Align short single-end reads.",
    "min_per_identity": "After the alignment step, mapped reads will be "
    "reclassified as unmapped if their identity percentage falls below this "
    "value. If not set, there is no reclassification.",
    "output_no_hits": "Report both matching and non-matching queries. "
    "WARNING: always use the default setting for this "
    "option unless if you know what you are doing! If "
    "you set this option to False, your sequences and "
    "feature table will need to be filtered to exclude "
    "unclassified sequences, otherwise you may run into "
    "errors downstream from missing feature IDs. Set to "
    "True to mirror default Minimap2 search.",
}
minimap2_search_params = {
    "n_threads": Int % Range(1, None),
    "maxaccepts": Int % Range(1, None),
    "preset": Str % Choices(["map-ont", "map-hifi", "map-pb", "sr"]),
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "output_no_hits": Bool,
}
minimap2_search_dsc = (
    "Search for top hits in a reference database using alignment between the "
    "query sequences and reference database sequences using Minimap2. Returns a "
    "report of the top M hits for each query (where M=maxaccepts)."
)

# classify-consensus-minimap2
classify_consensus_minimap2_inputs = {
    "query": FeatureData[Sequence],
    "index": Minimap2IndexDB,
    "reference": FeatureData[Sequence],
    "reference_taxonomy": FeatureData[Taxonomy],
}
classify_consensus_minimap2_outputs = [
    ("search_results", FeatureData[PairwiseAlignmentMN2]),
    ("classification", FeatureData[Taxonomy]),
]
classify_consensus_minimap2_inputs_dsc = {
    "query": "Query sequences.",
    "index": "Minimap2 indexed database. " "Incompatible with reference.",
    "reference": "Reference sequences. Incompatible with index.",
    "reference_taxonomy": "Reference taxonomy labels.",
}

classify_consensus_minimap2_outputs_dsc = {
    "search_results": "Top hits for each query.",
    "classification": "Taxonomy classifications of query sequences.",
}
classify_consensus_minimap2_params = {
    "maxaccepts": Int % Range(1, None),
    "preset": Str % Choices(["map-ont", "map-hifi", "map-pb", "sr"]),
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "output_no_hits": Bool,
    "n_threads": Int % Range(1, None),
    "min_consensus": Float % Range(0.5, 1.0, inclusive_end=True, inclusive_start=False),
    "unassignable_label": Str,
}
classify_consensus_minimap2_param_dsc = {
    "n_threads": "Number of threads to use.",
    "maxaccepts": (
        "Maximum number of hits to keep for each query. Minimap2 will "
        "choose the first N hits in the reference database that "
        "exceed perc_identity similarity to query."
    ),
    "preset": "The preset parameter applies multiple options at the same time "
    "during the mapping process of Minimap2. 1) map-ont: Align noisy long reads "
    "of ~10% error rate to a reference genome. 2) map-hifi: Align PacBio "
    "high-fidelity (HiFi) reads to a reference genome. 3) map-pb: Align "
    "older PacBio continuous long reads (CLR) to a reference genome. "
    "4) sr: Align short single-end reads.",
    "min_per_identity": "After the alignment step, mapped reads will be "
    "reclassified as unmapped if their identity percentage falls below this "
    "value. If not set, there is no reclassification.",
    "min_consensus": "Minimum fraction of assignments must match top "
    "hit to be accepted as consensus assignment.",
    "unassignable_label": "Annotation given to sequences without any hits.",
}
classify_consensus_minimap2_dsc = (
    "Assign taxonomy to query sequences using Minimap2. Performs "
    "alignment between query and reference reads, then "
    "assigns consensus taxonomy to each query sequence."
)

# find-consensus-annotation
find_consensus_annotation_inputs = {
    "search_results": FeatureData[PairwiseAlignmentMN2],
    "reference_taxonomy": FeatureData[Taxonomy],
}
find_consensus_annotation_params = {
    "min_consensus": Float % Range(0.5, 1.0, inclusive_end=True, inclusive_start=False),
    "unassignable_label": Str,
}
find_consensus_annotation_params_dsc = {
    "min_consensus": "Minimum fraction of assignments must match top "
    "hit to be accepted as consensus assignment.",
    "unassignable_label": "Annotation given when no consensus is found.",
}
find_consensus_annotation_outputs = [("consensus_taxonomy", FeatureData[Taxonomy])]
find_consensus_annotation_inputs_dsc = {
    "search_results": "Search results in PairwiseAlignmentMN2 output format",
    "reference_taxonomy": "Reference taxonomy labels.",
}
find_consensus_annotation_outputs_dsc = {
    "consensus_taxonomy": "Consensus taxonomy and scores."
}
find_consensus_annotation_dsc = (
    "Find consensus annotation for each query searched against "
    "a reference database, by finding the least common ancestor "
    "among one or more semicolon-delimited hierarchical "
    "annotations. Note that the annotation hierarchy is assumed "
    "to have an even number of ranks."
)
