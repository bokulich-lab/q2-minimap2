# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from q2_types.metadata import ImmutableMetadata
from q2_types.per_sample_sequences import (
    AlignmentMap,
    PairedEndSequencesWithQuality,
    SequencesWithQuality,
)
from q2_types.sample_data import SampleData
from qiime2.plugin import Bool, Choices, Float, Int, Range, Str, Threads, TypeMatch

from q2_minimap2.types._type import Minimap2IndexDB, PairwiseAlignmentMN2

T = TypeMatch([SequencesWithQuality, PairedEndSequencesWithQuality])

# Minimap2 presets suitable for searching reads against a reference database.
# The splice and ava presets are deliberately left out: ava sets -X, which skips
# the query-to-reference mappings every action here depends on, and splice
# introduces N CIGAR operations that the identity calculation does not count,
# which would let a spurious long-range alignment pass an identity threshold.
# asm5 and asm10 are omitted because they reject the divergence routinely seen
# between environmental reads and a reference database.
MAPPING_PRESETS = [
    "map-ont",
    "lr:hq",
    "map-hifi",
    "map-pb",
    "map-iclr",
    "asm20",
    "sr",
]

PRESET_DESCRIPTION = (
    "The preset parameter applies multiple options at the same time "
    "during the mapping process of Minimap2. 1) map-ont: Align noisy long reads "
    "of ~10% error rate to a reference genome. 2) lr:hq: Align accurate long "
    "reads of <1% error rate, such as Nanopore Q20+ data. 3) map-hifi: Align "
    "PacBio high-fidelity (HiFi) reads to a reference genome. 4) map-pb: Align "
    "older PacBio continuous long reads (CLR) to a reference genome. "
    "5) map-iclr: Align Illumina Complete Long Reads. 6) asm20: Align sequences "
    "diverged by up to roughly 5%, such as assembled contigs. 7) sr: Align "
    "short single-end reads. Note that lr:hq and map-hifi tolerate less "
    "divergence than map-ont, so map-ont remains the safer choice when the "
    "reference database may only contain distant relatives."
)

# filter_reads
filter_reads_inputs = {
    "query": SampleData[T],
    "index": Minimap2IndexDB,
    "reference": FeatureData[Sequence],
}
filter_reads_outputs = [
    ("filtered_query", SampleData[T]),
    ("filter_stats", ImmutableMetadata),
]
filter_reads_inputs_dsc = {
    "query": "Sequences to be filtered.",
    "index": "Minimap2 index database. Incompatible with reference.",
    "reference": "Reference sequences. Incompatible with index.",
}
filter_reads_outputs_dsc = {
    "filtered_query": "The resulting filtered sequences.",
    "filter_stats": "Per-sample counts of the reads read, kept and removed.",
}
filter_reads_params = {
    "n_threads": Threads,
    "preset": Str % Choices(MAPPING_PRESETS),
    "keep": Str % Choices(["mapped", "unmapped"]),
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "matching_score": Int,
    "mismatching_penalty": Int,
    "gap_open_penalty": Int % Range(1, None),
    "gap_extension_penalty": Int % Range(1, None),
}
filter_reads_param_dsc = {
    "n_threads": "Number of threads to use. Use 'auto' to use all available cores.",
    "preset": PRESET_DESCRIPTION,
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
    "n_threads": Threads,
    "preset": Str % Choices(MAPPING_PRESETS),
    "extract": Str % Choices(["mapped", "unmapped"]),
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "matching_score": Int,
    "mismatching_penalty": Int,
    "gap_open_penalty": Int % Range(1, None),
    "gap_extension_penalty": Int % Range(1, None),
}
extract_reads_param_dsc = {
    "n_threads": "Number of threads to use. Use 'auto' to use all available cores.",
    "preset": PRESET_DESCRIPTION,
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
    "preset": Str % Choices(MAPPING_PRESETS),
}
build_index_param_dsc = {
    "preset": "This option applies multiple settings at the same time during "
    "the indexing process. This value must match the mapping preset used in "
    "the actions that consume the index: Minimap2 keeps the k-mer and window "
    "settings baked into the index and silently ignores the ones implied by a "
    "different mapping preset, which changes the alignments that are found. "
    + PRESET_DESCRIPTION,
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
    "n_threads": "Number of threads to use. Use 'auto' to use all available cores.",
    "maxaccepts": "Maximum number of hits to keep for each query. When "
    "min_per_identity is set, the identity filter is applied first, so this "
    "keeps the top N of the hits that already satisfy it.",
    "preset": PRESET_DESCRIPTION,
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
    "n_threads": Threads,
    "maxaccepts": Int % Range(1, None),
    "preset": Str % Choices(MAPPING_PRESETS),
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
    "preset": Str % Choices(MAPPING_PRESETS),
    "min_per_identity": Float % Range(0.0, 1.0, inclusive_end=True),
    "output_no_hits": Bool,
    "n_threads": Threads,
    "min_consensus": Float % Range(0.5, 1.0, inclusive_end=True, inclusive_start=False),
    "unassignable_label": Str,
}
classify_consensus_minimap2_param_dsc = {
    "n_threads": "Number of threads to use. Use 'auto' to use all available cores.",
    "maxaccepts": (
        "Maximum number of hits to keep for each query. The identity filter "
        "is applied first, so this keeps the top N of the hits that already "
        "satisfy min_per_identity."
    ),
    "preset": PRESET_DESCRIPTION,
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


# align
align_inputs = {
    "query": SampleData[T],
    "index": Minimap2IndexDB,
    "reference": FeatureData[Sequence],
}
align_outputs = [("alignment", SampleData[AlignmentMap])]
align_inputs_dsc = {
    "query": "Sequences to align.",
    "index": "Minimap2 index database. Incompatible with reference.",
    "reference": "Reference sequences. Incompatible with index.",
}
align_outputs_dsc = {
    "alignment": "Coordinate-sorted BAM alignments, one per sample.",
}
align_params = {
    "n_threads": Threads,
    "preset": Str % Choices(MAPPING_PRESETS),
    "matching_score": Int,
    "mismatching_penalty": Int,
    "gap_open_penalty": Int % Range(1, None),
    "gap_extension_penalty": Int % Range(1, None),
}
align_param_dsc = {
    "n_threads": "Number of threads to use. Use 'auto' to use all " "available cores.",
    "preset": PRESET_DESCRIPTION,
    "matching_score": "Matching score.",
    "mismatching_penalty": "Mismatching penalty.",
    "gap_open_penalty": "Gap open penalty.",
    "gap_extension_penalty": "Gap extension penalty.",
}
align_dsc = (
    "Align sequencing reads to a set of reference sequences with Minimap2 and "
    "keep the alignment itself, rather than only the reads it selects. The "
    "result is a coordinate-sorted BAM per sample, which can be summarised "
    "with alignment-stats or handed to any other tool that reads BAM."
)

# alignment-stats
alignment_stats_inputs = {"alignments": SampleData[AlignmentMap]}
alignment_stats_outputs = [("stats", ImmutableMetadata)]
alignment_stats_inputs_dsc = {
    "alignments": "Alignments to summarise, as produced by align.",
}
alignment_stats_outputs_dsc = {
    "stats": "One row per sample, holding read counts, mapping rate, mean "
    "mapping quality and reference coverage.",
}
alignment_stats_dsc = (
    "Summarise Minimap2 alignments per sample. Reports how many reads mapped, "
    "how confidently they mapped, and how much of the reference they covered, "
    "as metadata that can be tabulated, plotted or joined to a sample "
    "metadata file."
)
