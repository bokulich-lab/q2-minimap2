# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
from unittest.mock import patch

import qiime2
from q2_types.feature_data import DNAFASTAFormat
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2.align import align
from q2_minimap2.alignment_stats import alignment_stats
from q2_minimap2.types._format import Minimap2IndexDBDirFmt

from .test_minimap2 import MinimapTestsBase

SAMPLE_IDS = {"sample_a", "sample_b", "sample_c"}


class TestAlign(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.single_end = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/single_end/"), mode="r"
        )
        self.paired_end = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/paired_end/"), mode="r"
        )
        self.index = Minimap2IndexDBDirFmt(
            self.get_data_path("filter_reads/index/"), mode="r"
        )
        self.reference = DNAFASTAFormat(
            self.get_data_path("filter_reads/dna-sequences.fasta"), mode="r"
        )

    def _bam_names(self, alignments):
        return {path.stem for path in alignments.path.glob("*.bam")}

    def test_align_single_end_writes_one_bam_per_sample(self):
        alignments = align(query=self.single_end, index=self.index)

        self.assertEqual(self._bam_names(alignments), SAMPLE_IDS)

    def test_align_paired_end_writes_one_bam_per_sample(self):
        alignments = align(query=self.paired_end, index=self.index)

        self.assertEqual(self._bam_names(alignments), SAMPLE_IDS)

    def test_align_accepts_a_reference_instead_of_an_index(self):
        alignments = align(query=self.single_end, reference=self.reference)

        self.assertEqual(self._bam_names(alignments), SAMPLE_IDS)

    def test_align_output_validates_as_bam(self):
        # BAMFormat runs `samtools quickcheck`, so this fails on a truncated or
        # otherwise malformed BAM
        alignments = align(query=self.single_end, index=self.index)
        alignments.validate()

    def test_align_rejects_both_reference_and_index(self):
        with self.assertRaisesRegex(ValueError, "Only one.*can be provided.*"):
            align(query=self.single_end, index=self.index, reference=self.reference)

    def test_align_requires_a_reference_or_an_index(self):
        with self.assertRaisesRegex(ValueError, "Either.*must be provided.*"):
            align(query=self.single_end)


class TestAlignmentStats(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.single_end = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/single_end/"), mode="r"
        )
        self.index = Minimap2IndexDBDirFmt(
            self.get_data_path("filter_reads/index/"), mode="r"
        )
        self.alignments = align(query=self.single_end, index=self.index)

    def test_one_row_per_sample(self):
        stats = alignment_stats(self.alignments)

        self.assertIsInstance(stats, qiime2.Metadata)
        self.assertEqual(set(stats.to_dataframe().index), SAMPLE_IDS)

    def test_reports_the_expected_columns(self):
        frame = alignment_stats(self.alignments).to_dataframe()

        for column in (
            "total_reads",
            "mapped_reads",
            "unmapped_reads",
            "percent_mapped",
            "mean_mapping_quality",
            "references_covered",
        ):
            self.assertIn(column, frame.columns)

    def test_counts_are_internally_consistent(self):
        frame = alignment_stats(self.alignments).to_dataframe()

        for _, row in frame.iterrows():
            self.assertEqual(
                row["mapped_reads"] + row["unmapped_reads"], row["total_reads"]
            )
            self.assertGreater(row["total_reads"], 0)
            self.assertTrue(0 <= row["percent_mapped"] <= 100)


class TestThreadResolution(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.single_end = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("filter_reads/single_end/"), mode="r"
        )
        self.index = Minimap2IndexDBDirFmt(
            self.get_data_path("filter_reads/index/"), mode="r"
        )

    def _thread_args(self, commands):
        return [cmd[cmd.index("-t") + 1] for cmd in commands if "-t" in cmd]

    @patch("q2_minimap2.align.run_cmd")
    @patch("q2_minimap2.align.get_available_cores", return_value=7)
    def test_zero_threads_resolves_to_every_available_core(self, cores, run_cmd):
        # The Threads type turns "auto" into 0 before the action sees it
        align(query=self.single_end, index=self.index, n_threads=0)

        cores.assert_called()
        commands = [call.args[0] for call in run_cmd.call_args_list]
        self.assertTrue(self._thread_args(commands))
        for value in self._thread_args(commands):
            self.assertEqual(value, "7")

    @patch("q2_minimap2.align.run_cmd")
    @patch("q2_minimap2.align.get_available_cores", return_value=7)
    def test_explicit_thread_count_is_passed_through(self, cores, run_cmd):
        align(query=self.single_end, index=self.index, n_threads=2)

        cores.assert_not_called()
        commands = [call.args[0] for call in run_cmd.call_args_list]
        for value in self._thread_args(commands):
            self.assertEqual(value, "2")


if __name__ == "__main__":
    unittest.main()
