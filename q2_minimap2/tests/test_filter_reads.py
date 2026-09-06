# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import os
import unittest

import qiime2
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2.filter_reads import filter_reads
from q2_minimap2.types._format import Minimap2IndexDBDirFmt

from .test_minimap2 import MinimapTestsBase

seq_ids_mapped = [
    "SARS2:6:73:941:1973#",
    "SARS2:6:73:231:3321#",
    "SARS2:6:73:233:3421#",
    "SARS2:6:73:552:2457#",
    "SARS2:6:73:567:7631#",
]
seq_ids_unmapped = ["SARS2:6:73:356:9806#"]

perc_id_mapped = [
    "SARS2:6:73:231:3321#",
    "SARS2:6:73:233:3421#",
    "SARS2:6:73:552:2457#",
    "SARS2:6:73:567:7631#",
]

perc_id_unmapped = ["SARS2:6:73:941:1973#", "SARS2:6:73:356:9806#"]


class TestFilterSingleEndReads(MinimapTestsBase):
    def setUp(self):
        super().setUp()

        self.query_single_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/single_end/"), mode="r"
        )
        self.query_paired_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/paired_end/"), mode="r"
        )
        self.minimap2_index = Minimap2IndexDBDirFmt(
            self.get_data_path("filter_reads/index/"), mode="r"
        )
        self.reference_reads = DNAFASTAFormat(
            self.get_data_path("filter_reads/dna-sequences.fasta"), mode="r"
        )

    def _check_ids(self, obs_seqs, included_ids, excluded_ids):

        fastq_files = [f for f in os.listdir(str(obs_seqs)) if f.endswith(".fastq.gz")]

        # Process each FASTQ.GZ file
        for obs_fp in fastq_files:
            file_path = os.path.join(str(obs_seqs), obs_fp)
            with gzip.open(file_path, "rt") as obs_fh:
                # Ensure the file is not empty
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                obs_fh.seek(0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    obs_seq_h, obs_seq, _, obs_qual = records
                    # Make sure seqs that map to genome were removed
                    obs_id = obs_seq_h.strip("@/012\n")
                    self.assertTrue(obs_id in included_ids)
                    self.assertTrue(obs_id not in excluded_ids)

    def test_filter_single_end_keep_unmapped(self):
        obs_seqs, _ = filter_reads(
            query=self.query_single_reads,
            index=self.minimap2_index,
            keep="unmapped",
        )

        self._check_ids(obs_seqs, seq_ids_unmapped, seq_ids_mapped)

    def test_filter_single_end_keep_mapped(self):
        obs_seqs, _ = filter_reads(
            query=self.query_single_reads,
            index=self.minimap2_index,
        )
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)

    def test_filter_single_end_keep_mapped_sr(self):
        obs_seqs, _ = filter_reads(
            query=self.query_single_reads,
            index=self.minimap2_index,
            preset="sr",
        )
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)

    def test_filter_single_end_keep_mapped_using_ref(self):
        obs_seqs, _ = filter_reads(
            query=self.query_single_reads,
            reference=self.reference_reads,
        )
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)

    def test_filter_single_end_keep_unmapped_with_perc_id(self):
        obs_seqs, _ = filter_reads(
            query=self.query_single_reads,
            index=self.minimap2_index,
            keep="unmapped",
            min_per_identity=0.99,
        )
        self._check_ids(obs_seqs, perc_id_unmapped, perc_id_mapped)

    def test_filter_single_end_keep_mapped_with_perc_id(self):
        obs_seqs, _ = filter_reads(
            query=self.query_single_reads,
            index=self.minimap2_index,
            keep="mapped",
            min_per_identity=0.99,
        )
        self._check_ids(obs_seqs, perc_id_mapped, perc_id_unmapped)

    def test_both_reference_and_index_provided(self):
        with self.assertRaises(ValueError) as context:
            filter_reads(
                query=self.query_single_reads,
                index=self.minimap2_index,
                reference=self.reference_reads,
            )
        self.assertIn(
            "Only one of reference or index can be provided",
            str(context.exception),
        )

    def test_neither_reference_nor_index_provided(self):
        with self.assertRaises(ValueError) as context:
            filter_reads(
                query=self.query_single_reads,
                index=None,
                reference=None,
            )
        self.assertIn(
            "Either reference or index must be provided",
            str(context.exception),
        )

    def test_filter_paired_end_keep_unmapped(self):
        obs_seqs, _ = filter_reads(
            query=self.query_paired_reads,
            index=self.minimap2_index,
            keep="unmapped",
        )
        self._check_ids(obs_seqs, seq_ids_unmapped, seq_ids_mapped)

    def test_filter_paired_end_keep_mapped(self):
        obs_seqs, _ = filter_reads(
            query=self.query_paired_reads,
            index=self.minimap2_index,
        )
        self._check_ids(obs_seqs, seq_ids_mapped, seq_ids_unmapped)


class TestFilterStats(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.query_single_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/single_end/"), mode="r"
        )
        self.minimap2_index = Minimap2IndexDBDirFmt(
            self.get_data_path("filter_reads/index/"), mode="r"
        )

    def test_stats_reported_for_every_sample(self):
        _, stats = filter_reads(
            query=self.query_single_reads, index=self.minimap2_index
        )

        self.assertIsInstance(stats, qiime2.Metadata)
        self.assertEqual(
            set(stats.to_dataframe().index), {"sample_a", "sample_b", "sample_c"}
        )

    def test_kept_and_removed_account_for_the_input(self):
        _, stats = filter_reads(
            query=self.query_single_reads, index=self.minimap2_index
        )

        for _, row in stats.to_dataframe().iterrows():
            self.assertEqual(
                row["retained_reads"] + row["removed_reads"], row["input_reads"]
            )
            self.assertGreater(row["input_reads"], 0)

    def test_paired_end_stats_match_the_reads_actually_written(self):
        # The paired path drops a mate whose partner was filtered out, so the
        # count must come from what was written, not from what passed the filter
        paired = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/paired_end/"), mode="r"
        )
        obs_seqs, stats = filter_reads(
            query=paired, index=self.minimap2_index, min_per_identity=0.99
        )

        frame = stats.to_dataframe()
        for sample_id, row in frame.iterrows():
            forward = os.path.join(
                str(obs_seqs), f"{sample_id}_S01_L001_R1_001.fastq.gz"
            )
            matching = [
                f
                for f in os.listdir(str(obs_seqs))
                if f.startswith(sample_id) and "R1" in f
            ]
            if not matching:
                continue
            forward = os.path.join(str(obs_seqs), matching[0])
            with gzip.open(forward, "rt") as handle:
                forward_reads = sum(1 for _ in handle) // 4
            reverse = forward.replace("R1", "R2")
            with gzip.open(reverse, "rt") as handle:
                reverse_reads = sum(1 for _ in handle) // 4

            self.assertEqual(row["retained_reads"], forward_reads + reverse_reads)

    def test_keeping_mapped_and_unmapped_partitions_the_input(self):
        # Every read is either kept as mapped or kept as unmapped, so the two
        # runs must add up to the input
        _, mapped = filter_reads(
            query=self.query_single_reads, index=self.minimap2_index, keep="mapped"
        )
        _, unmapped = filter_reads(
            query=self.query_single_reads, index=self.minimap2_index, keep="unmapped"
        )

        mapped_frame = mapped.to_dataframe()
        unmapped_frame = unmapped.to_dataframe()
        for sample_id, row in mapped_frame.iterrows():
            self.assertEqual(
                row["retained_reads"] + unmapped_frame.loc[sample_id]["retained_reads"],
                row["input_reads"],
            )


if __name__ == "__main__":
    unittest.main()
