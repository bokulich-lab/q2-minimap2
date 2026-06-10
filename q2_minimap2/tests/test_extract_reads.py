# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import unittest

from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2.extract_reads import extract_reads
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


class TestExtractreads(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.query1_reads = DNAFASTAFormat(
            self.get_data_path("extract_reads/extract_reads_test_input.fasta"), mode="r"
        )
        self.query2_reads = DNAFASTAFormat(
            self.get_data_path("extract_reads/query-reads.fasta"), mode="r"
        )
        self.minimap2_index = Minimap2IndexDBDirFmt(
            self.get_data_path("extract_reads/index/"), mode="r"
        )
        self.ref = DNAFASTAFormat(
            self.get_data_path("extract_reads/se-dna-sequences.fasta"), mode="r"
        )

    def _check_extracted_ids(self, extracted_sequences, included_ids, excluded_ids):
        with open(str(extracted_sequences), "rt") as obs_fh:
            self.assertNotEqual(len(obs_fh.readlines()), 0)
            obs_fh.seek(0)
            for records in itertools.zip_longest(*[obs_fh] * 4):
                obs_seq_h, obs_seq, _, obs_qual = records
                obs_id = obs_seq_h.strip(">/012\n")
                self.assertTrue(obs_id in included_ids)
                self.assertTrue(obs_id not in excluded_ids)

    def test_extract_unmapped(self):
        extracted_sequences = extract_reads(
            self.query1_reads, self.minimap2_index, extract="unmapped"
        )
        self._check_extracted_ids(extracted_sequences, seq_ids_unmapped, seq_ids_mapped)

    def test_extract_mapped(self):
        extracted_sequences = extract_reads(self.query1_reads, self.minimap2_index)
        self._check_extracted_ids(extracted_sequences, seq_ids_mapped, seq_ids_unmapped)

    def test_extract_unmapped_with_perc_id(self):
        extracted_sequences = extract_reads(
            self.query1_reads,
            self.minimap2_index,
            extract="unmapped",
            min_per_identity=0.99,
        )
        self._check_extracted_ids(extracted_sequences, perc_id_unmapped, perc_id_mapped)

    def test_extract_mapped_with_perc_id(self):
        extracted_sequences = extract_reads(
            self.query1_reads,
            self.minimap2_index,
            extract="mapped",
            min_per_identity=0.99,
        )
        self._check_extracted_ids(extracted_sequences, perc_id_mapped, perc_id_unmapped)

    def test_extract_mapped_using_reference(self):
        extracted_sequences = extract_reads(self.query2_reads, reference=self.ref)
        obs_fp = str(extracted_sequences)
        correct_output_fp = self.get_data_path("extract_reads/extracted_mapped.fasta")
        with open(correct_output_fp, "r") as file1, open(obs_fp, "r") as file2:
            true_fasta_content = file1.read()
            output_fasta_content = file2.read()
            self.assertEqual(true_fasta_content, output_fasta_content)

    def test_extract_reads_both_ref_and_index(self):
        with self.assertRaisesRegex(ValueError, "Only one.*can be provided.*"):
            extract_reads(
                self.query2_reads,
                reference=self.ref,
                index=self.minimap2_index,
            )
        with self.assertRaisesRegex(ValueError, "Either.*must be provided.*"):
            extract_reads(self.query2_reads)


if __name__ == "__main__":
    unittest.main()
