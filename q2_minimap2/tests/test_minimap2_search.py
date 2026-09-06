# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest

import pandas as pd
from pandas.testing import assert_frame_equal
from q2_types.per_sample_sequences import CasavaOneEightSingleLanePerSampleDirFmt

from q2_minimap2.minimap2_search import (
    construct_command,
    filter_by_maxaccepts,
    filter_by_perc_identity,
    minimap2_search,
)
from q2_minimap2.types._format import Minimap2IndexDBDirFmt, read_paf

from .test_minimap2 import MinimapTestsBase


class TestFilterByMaxAccepts(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        # Retrieve the paths to test data files
        self.initial_PairwiseAlignmentMN2_file = self.get_data_path(
            "minimap2_search/initial_paf_file.paf"
        )
        self.exp_PairwiseAlignmentMN2_f_max1 = self.get_data_path(
            "minimap2_search/expected_paf_file_max1.paf"
        )
        self.exp_PairwiseAlignmentMN2_f_max2 = self.get_data_path(
            "minimap2_search/expected_paf_file_max2.paf"
        )
        self.exp_PairwiseAlignmentMN2_f_max3 = self.get_data_path(
            "minimap2_search/expected_paf_file_max3.paf"
        )

    def test_filter_by_maxaccepts_max1(self):
        self._test_filter_by_maxaccepts(1, self.exp_PairwiseAlignmentMN2_f_max1)

    def test_filter_by_maxaccepts_max2(self):
        self._test_filter_by_maxaccepts(2, self.exp_PairwiseAlignmentMN2_f_max2)

    def test_filter_by_maxaccepts_max3(self):
        self._test_filter_by_maxaccepts(3, self.exp_PairwiseAlignmentMN2_f_max3)

    def _test_filter_by_maxaccepts(self, max_accepts, expected_file_path):
        # Load input and expected data
        input_df = pd.read_csv(
            self.initial_PairwiseAlignmentMN2_file, sep="\t", header=None
        )
        expected_df = pd.read_csv(expected_file_path, sep="\t", header=None)

        # Generate results from function
        result_df = filter_by_maxaccepts(input_df, max_accepts)
        result_df.reset_index(drop=True, inplace=True)

        # Assert that the two data frames are equal
        assert_frame_equal(result_df, expected_df)


class TestFilterByPercIdentity(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.PairwiseAlignmentMN2_file = self.get_data_path(
            "minimap2_search/initial_paf_file.paf"
        )
        self.expected_PairwiseAlignmentMN2_file_perc_85 = self.get_data_path(
            "minimap2_search/expected_paf_file_perc_85.paf"
        )
        self.expected_PairwiseAlignmentMN2_file_perc_80 = self.get_data_path(
            "minimap2_search/expected_paf_file_perc_80.paf"
        )

    def test_filter_by_perc_identity_85(self):
        self._test_filter_by_perc_identity(
            0.85, self.expected_PairwiseAlignmentMN2_file_perc_85
        )

    def test_filter_by_perc_identity_80(self):
        self._test_filter_by_perc_identity(
            0.8, self.expected_PairwiseAlignmentMN2_file_perc_80
        )

    def _test_filter_by_perc_identity(self, perc_identity, expected_file_path):
        # Load input and expected data
        input_df = pd.read_csv(self.PairwiseAlignmentMN2_file, sep="\t", header=None)
        expected_df = pd.read_csv(expected_file_path, sep="\t", header=None)

        # Generate results from function
        result_df = filter_by_perc_identity(input_df, perc_identity, True)
        result_df.reset_index(drop=True, inplace=True)

        # Assert that the two data frames are equal
        assert_frame_equal(result_df, expected_df)


class TestMinimap2(MinimapTestsBase):
    def setUp(self):
        super().setUp()

        self.query_reads = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("minimap2_search/query_seqs.fasta"), mode="r"
        )
        self.index_database = Minimap2IndexDBDirFmt(
            self.get_data_path("minimap2_search/minimap2_test_index/"), mode="r"
        )
        self.ref = CasavaOneEightSingleLanePerSampleDirFmt(
            self.get_data_path("minimap2_search/se-dna-sequences.fasta"), mode="r"
        )

    def test_minimap2(self):
        # Perform the minimap2 search and store the result in a DataFrame
        search_results_df = minimap2_search(self.query_reads, self.index_database)

        # Define the path to the expected output file
        expected_output_path = self.get_data_path(
            "minimap2_search/minimap2_test_paf.paf"
        )

        # Read the expected output content from the file
        with open(expected_output_path, "r") as file:
            expected_output_content = file.read().strip()

        # Convert the DataFrame to a tab-separated string without index and header
        results_content = search_results_df.to_csv(
            sep="\t", index=False, header=False
        ).strip()
        # Clean the extra tabs
        results_content_cleaned = "\n".join(
            [
                "\t".join(line.split()).rstrip("\t")
                for line in results_content.split("\n")
            ]
        )

        # Assert that the search results match the expected output
        self.assertEqual(expected_output_content, results_content_cleaned)

    def test_minimap2_only_hits(self):
        search_results_df = minimap2_search(
            self.query_reads, self.index_database, output_no_hits=False
        )
        expected_output_path = self.get_data_path(
            "minimap2_search/minimap2_only_hits_test_paf.paf"
        )

        with open(expected_output_path, "r") as file:
            expected_output_content = file.read()
            results_content = search_results_df.to_csv(
                sep="\t", index=False, header=False
            )
            # Assert that the search results match the expected output
            self.assertEqual(expected_output_content, results_content)

    def test_minimap2_output_consistency(self):
        result1 = minimap2_search(
            self.query_reads, self.index_database, output_no_hits=False
        )
        result2 = minimap2_search(
            self.query_reads, reference=self.ref, output_no_hits=False
        )
        self.assertEqual(
            result1.to_csv(sep="\t", index=False, header=False),
            result2.to_csv(sep="\t", index=False, header=False),
        )

    def test_minimap2_both_ref_and_index(self):
        with self.assertRaisesRegex(ValueError, "Only one.*can be provided.*"):
            minimap2_search(
                self.query_reads,
                reference=self.ref,
                index=self.index_database,
            )
        with self.assertRaisesRegex(ValueError, "Either.*must be provided.*"):
            minimap2_search(self.query_reads)


class TestFilterByPercIdentityNoHits(MinimapTestsBase):
    # Builds a PAF frame where every query is the same length, which is what a
    # fixed-length amplicon run looks like
    def _fixed_length_paf(self, n_queries, matches=100):
        rows = [
            [f"query{i}", 150, 0, 150, "+", f"ref{i}", 1500, 0, 150, matches, 150, 60]
            for i in range(n_queries)
        ]
        return pd.DataFrame(rows)

    def test_no_hit_rows_kept_per_query_not_per_length(self):
        # Every query fails the threshold, so every one of them should come
        # back as its own no-hit row rather than collapsing into a single row
        df = self._fixed_length_paf(5)
        result = filter_by_perc_identity(df, 0.95, True)

        self.assertEqual(len(result), 5)
        self.assertEqual(sorted(result[0]), [f"query{i}" for i in range(5)])

    def test_no_hit_row_not_added_for_query_with_accepted_hit(self):
        # One query has a passing and a failing alignment; only the passing one
        # should survive, with no fabricated unmapped row alongside it
        df = pd.DataFrame(
            [
                ["q1", 150, 0, 150, "+", "refA", 1500, 0, 150, 149, 150, 60],
                ["q1", 150, 0, 150, "+", "refB", 1500, 0, 150, 100, 150, 60],
            ]
        )
        result = filter_by_perc_identity(df, 0.95, True)

        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0][5], "refA")

    def test_numeric_reference_ids_accept_no_hit_marker(self):
        # Numeric reference identifiers must not stop the "*" marker from being
        # written into the strand and target columns
        df = pd.DataFrame(
            [["q1", 150, 0, 150, "+", "1111561", 1500, 0, 150, 100, 150, 60]]
        )
        result = filter_by_perc_identity(df, 0.95, True)

        self.assertEqual(result.iloc[0][4], "*")
        self.assertEqual(result.iloc[0][5], "*")


class TestReadPaf(MinimapTestsBase):
    def _write(self, lines):
        path = os.path.join(self.temp_dir.name, "test.paf")
        with open(path, "w") as file:
            file.write("".join(line + "\n" for line in lines))
        return path

    def test_rows_of_differing_width(self):
        # A no-hit record is narrower than an aligned one, and the narrow row
        # coming first must not fix the width for the whole file
        narrow = "\t".join(
            ["junk", "100", "0", "0", "*", "*", "0", "0", "0", "0", "0", "0"]
        )
        wide = "\t".join(
            ["q1", "75", "0", "75", "+", "r1", "150", "0", "75", "75", "75", "60"]
            + ["NM:i:0", "tp:A:P"]
        )
        df = read_paf(self._write([narrow, wide]))

        self.assertEqual(len(df), 2)
        self.assertEqual(len(df.columns), 14)

    def test_empty_file(self):
        # A run that produced no alignments at all leaves an empty file
        df = read_paf(self._write([]))

        self.assertEqual(len(df), 0)
        self.assertEqual(len(df.columns), 12)

    def test_identifiers_are_read_verbatim(self):
        # "NA" must stay a name rather than becoming a missing value, and a
        # zero-padded identifier must not be turned into a number
        row = "\t".join(
            ["NA", "75", "0", "75", "+", "007", "150", "0", "75", "75", "75", "60"]
        )
        df = read_paf(self._write([row]))

        self.assertEqual(df.iloc[0][0], "NA")
        self.assertEqual(df.iloc[0][5], "007")


class TestConstructCommand(MinimapTestsBase):
    def test_secondary_alignment_limit_follows_maxaccepts(self):
        # Minimap2 keeps only 5 secondary alignments unless told otherwise,
        # which would cap every query below a larger maxaccepts
        cmd = construct_command(
            "idx", "query.fasta", 1, "map-ont", "out.paf", False, 10
        )

        self.assertIn("-N", cmd)
        self.assertEqual(cmd[cmd.index("-N") + 1], "10")

    def test_split_prefix_is_passed(self):
        cmd = construct_command(
            "idx", "query.fasta", 1, "map-ont", "out.paf", False, 1, "/tmp/split"
        )

        self.assertIn("--split-prefix", cmd)
        self.assertEqual(cmd[cmd.index("--split-prefix") + 1], "/tmp/split")


if __name__ == "__main__":
    unittest.main()
