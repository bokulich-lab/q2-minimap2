# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from q2_minimap2._filtering_utils import (
    calculate_identity,
    collate_sam_inplace,
    convert_to_fasta,
    convert_to_fastq,
    get_alignment_length,
    make_mn2_cmd,
    process_paired_sam_flags,
    process_sam_file,
    run_cmd,
    set_penalties,
)

from .test_minimap2 import MinimapTestsBase


class TestSetPenalties(MinimapTestsBase):

    # Test with all parameters provided to ensure correct option string formation
    def test_all_arguments_provided(self):
        self.assertEqual(
            set_penalties(1, -1, 4, 1), ["-A", "1", "-B", "-1", "-O", "4", "-E", "1"]
        )

    # Test with some parameters as None to check function's handling of optional
    # arguments
    def test_some_arguments_none(self):
        # Expecting options only for mismatch and gap_e since others are None
        self.assertEqual(set_penalties(None, -1, None, 2), ["-B", "-1", "-E", "2"])
        # Expecting options only for match and gap_o since others are None
        self.assertEqual(set_penalties(1, None, 4, None), ["-A", "1", "-O", "4"])

    # Test with all parameters as None to verify it returns an empty list, indicating
    # no options
    def test_all_arguments_none(self):
        self.assertEqual(set_penalties(None, None, None, None), [])


class TestCalculateIdentity(MinimapTestsBase):
    # Test case with NM tag present
    def test_with_nm_tag(self):
        aln = "SOME\tDATA\tNM:i:5"
        total_length = 100
        expected_identity = 0.95  # (100 - 5) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Test case without NM tag, should default to 0 mismatches
    def test_without_nm_tag(self):
        aln = "SOME\tDATA"
        total_length = 100
        expected_identity = 1.0  # (100 - 0) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Test case with NM tag present but 0 mismatches
    def test_with_nm_tag_zero_mismatches(self):
        aln = "SOME\tDATA\tNM:i:0"
        total_length = 100
        expected_identity = 1.0  # (100 - 0) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Test case with total length 0 to check for division by zero handling
    def test_total_length_zero(self):
        aln = "SOME\tDATA\tNM:i:5"
        total_length = 0
        with self.assertRaises(ZeroDivisionError):
            calculate_identity(aln, total_length)

    # Test case with negative mismatches, which shouldn't happen but let's test
    def test_negative_mismatches(self):
        aln = "SOME\tDATA\tNM:i:-5"
        total_length = 100
        expected_identity = 1.05  # (100 - (-5)) / 100, though practically impossible
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Tests the function's response to a non-integer mismatch count in the NM tag
    def test_malformed_nm_tag(self):
        aln = "SOME\tDATA\tNM:i:notanumber"
        total_length = 100
        # Expecting a ValueError due to the inability to parse "notanumber" as
        # an integer
        with self.assertRaises(ValueError):
            calculate_identity(aln, total_length)

    # Verifies the function's behavior when multiple NM tags are present in
    # the alignment
    def test_multiple_nm_tags(self):
        aln = "SOME\tDATA\tNM:i:5\tANOTHER\tTAG\tNM:i:10"
        total_length = 100
        # Assuming the function picks the first NM tag
        expected_identity = 0.95  # (100 - 5) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Checks how the function calculates identity with mismatches exceeding the
    # total length
    def test_large_number_of_mismatches(self):
        aln = "SOME\tDATA\tNM:i:150"
        total_length = 100
        # More mismatches than total length, resulting in negative identity
        expected_identity = -0.5  # (100 - 150) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Assesses the function's ability to handle special characters in the
    # alignment string
    def test_special_characters_in_alignment(self):
        aln = "SOME\t$%^&*()\tDATA\tNM:i:5"
        total_length = 100
        expected_identity = 0.95  # (100 - 5) / 100
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)

    # Evaluates the function's handling of a floating-point number for the total
    # alignment length
    def test_floating_point_total_length(self):
        aln = "SOME\tDATA\tNM:i:5"
        total_length = 100.5
        # Test handling of floating point total lengths
        expected_identity = (100.5 - 5) / 100.5
        self.assertAlmostEqual(calculate_identity(aln, total_length), expected_identity)


class TestGetAlignmentLength(MinimapTestsBase):
    # Test with the CIGAR string indicating no alignment
    def test_no_alignment(self):
        self.assertEqual(get_alignment_length("*"), 0)

    # Test with a CIGAR string that has only matches
    def test_only_matches(self):
        self.assertEqual(get_alignment_length("10M"), 10)

    # Test with a CIGAR string that has only insertions
    def test_only_insertions(self):
        self.assertEqual(get_alignment_length("5I"), 5)

    # Test with a CIGAR string that has only deletions
    def test_only_deletions(self):
        self.assertEqual(get_alignment_length("3D"), 3)

    # Test with a CIGAR string that combines matches, insertions, and deletions
    def test_combination(self):
        self.assertEqual(get_alignment_length("5M2I3D"), 10)

    # Test with a more complex CIGAR string
    def test_complex_cigar(self):
        self.assertEqual(get_alignment_length("6M1I4M1D5M"), 17)

    # Test with large numbers in the CIGAR string
    def test_large_numbers(self):
        self.assertEqual(get_alignment_length("100M50I25D"), 175)

    # Test with a CIGAR string that includes operations not contributing to alignment
    #  length (e.g., 'S' for soft clipping)
    def test_with_non_contributing_operations(self):
        self.assertEqual(get_alignment_length("5M3S2I4D"), 11)

    # Test with an empty CIGAR string
    def test_empty_cigar(self):
        self.assertEqual(get_alignment_length(""), 0)

    # Test with invalid characters in the CIGAR string
    def test_invalid_characters(self):
        self.assertEqual(
            get_alignment_length("10M2X3M"), 13
        )  # Assuming 'X' is treated as invalid and ignored

    # Test CIGAR strings with missing numbers before operation characters
    def test_missing_numbers(self):
        self.assertEqual(
            get_alignment_length("MID"), 0
        )  # Assuming it results in 0 as no valid operations are detected

    # Test handling of leading and trailing non-operation characters
    def test_leading_trailing_non_operations(self):
        self.assertEqual(
            get_alignment_length("A10M3DZ"), 13
        )  # Assuming non-operations are ignored

    # Test sequential operations without numbers in-between
    def test_sequential_operations_without_numbers(self):
        self.assertEqual(get_alignment_length("2M3I2D2M"), 9)

    # Test operations specified with a length of 0
    def test_zero_length_operations(self):
        # Assuming 0-length operations are valid and result in 0 contribution
        self.assertEqual(get_alignment_length("0M10I0D"), 10)

    # Test with non-standard operations that don't contribute to alignment length
    def test_non_standard_operations(self):
        self.assertEqual(
            get_alignment_length("10M3S2H"), 10
        )  # Assuming 'S' and 'H' are ignored for alignment length


class TestProcessSamFile(MinimapTestsBase):
    def setUp(self):
        super().setUp()

        self.samfile = self.get_data_path("filtering_utils/initial.sam")
        self.only_primary = self.get_data_path("filtering_utils/mapped.sam")
        self.only_unmapped = self.get_data_path("filtering_utils/unmapped.sam")

    def test_keep_primary_mapped(self):
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            # Concatenate the temporary directory path and the original file name
            tmpfile = os.path.join(tmpdir, "initial_tmp.sam")

            # Copy the SAM file to the temporary location
            shutil.copy(self.samfile, tmpfile)

            process_sam_file(tmpfile, "mapped", None)

            # Check each line for equality using splitlines()
            with open(tmpfile, "r") as tmpfile_content, open(
                self.only_primary, "r"
            ) as primary_content:
                tmp_lines = tmpfile_content.read().splitlines()
                primary_lines = primary_content.read().splitlines()

                tmp_index = 0
                primary_index = 0

                while tmp_index < len(tmp_lines) and primary_index < len(primary_lines):
                    tmp_line = tmp_lines[tmp_index]
                    primary_line = primary_lines[primary_index]

                    if tmp_line.startswith("@PG"):
                        # Skip the line starting with '@PG' in tmp_lines
                        tmp_index += 1
                        primary_index += 1
                        continue

                    # Explicitly print the line
                    self.assertEqual(tmp_line, primary_line)

                    tmp_index += 1
                    primary_index += 1

    def test_keep_primary_mapped_with_perc_id(self):
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            # Concatenate the temporary directory path and the original file name
            tmpfile = os.path.join(tmpdir, "initial_tmp.sam")

            # Copy the SAM file to the temporary location
            shutil.copy(self.samfile, tmpfile)

            process_sam_file(tmpfile, "mapped", 0.3)

            # Check each line for equality using splitlines()
            with open(tmpfile, "r") as tmpfile_content, open(
                self.only_primary, "r"
            ) as primary_content:
                tmp_lines = tmpfile_content.read().splitlines()
                primary_lines = primary_content.read().splitlines()

                tmp_index = 0
                primary_index = 0

                while tmp_index < len(tmp_lines) and primary_index < len(primary_lines):
                    tmp_line = tmp_lines[tmp_index]
                    primary_line = primary_lines[primary_index]

                    if tmp_line.startswith("@PG"):
                        # Skip the line starting with '@PG' in tmp_lines
                        tmp_index += 1
                        primary_index += 1
                        continue

                    # Explicitly print the line
                    self.assertEqual(tmp_line, primary_line)

                    tmp_index += 1
                    primary_index += 1

    def test_keep_unmapped(self):
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            # Concatenate the temporary directory path and the original file name
            tmpfile = os.path.join(tmpdir, "initial_tmp.sam")

            # Copy the SAM file to the temporary location
            shutil.copy(self.samfile, tmpfile)

            process_sam_file(tmpfile, "unmapped", None)

            # Check each line for equality using splitlines()
            with open(tmpfile, "r") as tmpfile_content, open(
                self.only_unmapped, "r"
            ) as unmapped_content:
                tmp_lines = tmpfile_content.read().splitlines()
                unmapped_lines = unmapped_content.read().splitlines()

                tmp_index = 0
                unmapped_index = 0

                while tmp_index < len(tmp_lines) and unmapped_index < len(
                    unmapped_lines
                ):
                    tmp_line = tmp_lines[tmp_index]
                    unmapped_line = unmapped_lines[unmapped_index]

                    if tmp_line.startswith("@PG"):
                        # Skip the line starting with '@PG' in tmp_lines
                        tmp_index += 1
                        unmapped_index += 1
                        continue

                    # Explicitly print the line
                    self.assertEqual(tmp_line, unmapped_line)

                    tmp_index += 1
                    unmapped_index += 1


class TestConvertToFasta(MinimapTestsBase):
    def test_convert_to_fasta(self):
        # Mock input parameters
        _reads = "read1.sam"
        n_threads = 4
        samfile_filepath = "output.sam"

        # Expected command based on the provided inputs
        expected_cmd = [
            "samtools",
            "fasta",
            "-0",
            str(_reads),
            "-s",
            "/dev/null",
            "-@",
            str(n_threads),
            "-n",
            str(samfile_filepath),
        ]

        # Call the function and compare the result with the expected command
        actual_cmd = convert_to_fasta(_reads, n_threads, samfile_filepath)
        self.assertEqual(actual_cmd, expected_cmd)


class TestConvertToFastq(MinimapTestsBase):
    def test_convert_to_fasta(self):
        # Mock input parameters
        _reads = "read1.sam"
        n_threads = 4
        samfile_filepath = "output.sam"

        # Expected command based on the provided inputs
        expected_cmd = [
            "samtools",
            "fastq",
            *_reads,
            "-s",
            "/dev/null",
            "-@",
            str(n_threads),
            "-n",
            str(samfile_filepath),
        ]

        # Call the function and compare the result with the expected command
        actual_cmd = convert_to_fastq(_reads, n_threads, samfile_filepath, "single")
        self.assertEqual(actual_cmd, expected_cmd)


class TestMakeMn2Cmd(MinimapTestsBase):
    def test_make_mn2_cmd(self):
        # Mock input parameters
        mapping_preset = "map-preset"
        index = "/path/to/index"
        n_threads = 4
        penalties = ["-A", "1", "-B", "-1", "-O", "4", "-E", "1"]
        reads = "input.fastq"
        samf_fp = "output.sam"

        # Expected command based on the provided inputs
        expected_cmd = [
            "minimap2",
            "-a",
            "-x",
            mapping_preset,
            "/path/to/index",
            "-t",
            str(n_threads),
            "-o",
            "output.sam",
            "-A",
            "1",
            "-B",
            "-1",
            "-O",
            "4",
            "-E",
            "1",
            "input.fastq",
            "--split-prefix",
            "/tmp/split",
        ]

        # Call the function and compare the result with the expected command
        actual_cmd = make_mn2_cmd(
            mapping_preset,
            index,
            n_threads,
            penalties,
            reads,
            None,
            samf_fp,
            "/tmp/split",
        )
        self.assertEqual(actual_cmd, expected_cmd)


class TestConvertToFastqPaired(MinimapTestsBase):
    def test_convert_to_fastq_paired(self):
        _reads = ["read1.qc", "read2.fq"]
        n_threads = 4
        samfile_filepath = "output.sam"

        expected_cmd = [
            "samtools",
            "fastq",
            *_reads,
            "-0",
            "/dev/null",
            "-s",
            "/dev/null",
            "-@",
            str(n_threads),
            "-n",
            str(samfile_filepath),
        ]

        result_cmd = convert_to_fastq(_reads, n_threads, samfile_filepath, "paired")
        self.assertEqual(
            result_cmd,
            expected_cmd,
            "The generated command for converting SAM to FastQ paired-end does not "
            "match the expected command.",
        )


class TestMakeMn2PairedEndCmd(MinimapTestsBase):
    def test_make_mn2_cmd_paired(self):
        mapping_preset = "map-preset"
        index = "/path/to/index"
        n_threads = 4
        penalties = ["-A", "1", "-B", "-1", "-O", "4", "-E", "1"]
        reads1 = "input1.fastq"
        reads2 = "input2.fastq"
        samf_fp = "output.sam"

        expected_cmd = [
            "minimap2",
            "-a",
            "-x",
            mapping_preset,
            str(index),
            "-t",
            str(n_threads),
            "-o",
            samf_fp,
            "-A",
            "1",
            "-B",
            "-1",
            "-O",
            "4",
            "-E",
            "1",
            reads1,
            reads2,
        ]

        # Paired input deliberately omits --split-prefix: it would switch
        # Minimap2 into fragment mode and change which reads align
        result_cmd = make_mn2_cmd(
            mapping_preset,
            index,
            n_threads,
            penalties,
            reads1,
            reads2,
            samf_fp,
            "/tmp/split",
        )
        self.assertEqual(result_cmd, expected_cmd)


class TestRunCmd(MinimapTestsBase):
    @patch("subprocess.run")
    def test_run_cmd_success(self, mock_run):
        cmd = "samtools fastq -o output.fastq input.sam"
        description = "samtools command"
        mock_run.return_value = MagicMock(returncode=0)  # Simulate successful run

        # Attempt to run the command
        try:
            run_cmd(cmd, description)
        except Exception as e:
            self.fail(f"run_cmd raised an exception unexpectedly: {e}")

    @patch("subprocess.run")
    def test_run_cmd_exception(self, mock_run):
        """Test that run_cmd raises the correct exception when subprocess.run fails."""
        cmd = "samtools fastq -o output.fastq input.sam"
        description = "samtools command"
        mock_run.side_effect = subprocess.CalledProcessError(
            1, cmd
        )  # Simulate an error

        with self.assertRaises(Exception) as context:
            run_cmd(cmd, description)

        expected_message = (
            "An error was encountered while using samtools command, "
            "(return code 1), please inspect stdout and stderr to learn more."
        )
        self.assertEqual(str(context.exception), expected_message)


class TestCollateSamInplace(MinimapTestsBase):
    @patch("q2_minimap2._filtering_utils.run_cmd")
    @patch("shutil.move")
    def test_collate_sam_inplace(self, mock_shutil_move, mock_run_cmd):
        input_sam_path = "test_input.sam"
        temp_prefix = os.path.splitext(input_sam_path)[0] + "_temp_collate"
        output_sam_path = os.path.splitext(input_sam_path)[0] + "_collated.sam"

        collate_sam_inplace(input_sam_path)

        # Verify run_cmd was called with the correct command
        expected_cmd = [
            "samtools",
            "collate",
            "-u",
            "-o",
            output_sam_path,
            "-T",
            temp_prefix,
            input_sam_path,
        ]
        mock_run_cmd.assert_called_once_with(expected_cmd, "Samtools collate")

        # Verify shutil.move was called with the correct arguments
        mock_shutil_move.assert_called_once_with(output_sam_path, input_sam_path)


class TestProcessPairedSamFlags(MinimapTestsBase):
    def setUp(self):
        # Create a temporary SAM file for testing
        self.temp_dir = tempfile.mkdtemp()
        self.sam_file = os.path.join(self.temp_dir, "test.sam")
        self.create_sam_file()

    def tearDown(self):
        # Remove the temporary directory and file
        shutil.rmtree(self.temp_dir)

    def create_sam_file(self):
        # Mates of a pair share a read name, which is what the collated SAM the
        # function consumes actually looks like
        content = (
            "@HD\tVN:1.0\tSO:unsorted\n"
            "@SQ\tSN:chr1\tLN:248956422\n"
            "readA\t4\tchr1\t100\t255\t50M\t*\t0\t0\t*\t*\n"
            "readA\t4\tchr1\t150\t255\t50M\t*\t0\t0\t*\t*\n"
            "readB\t0\tchr1\t200\t255\t50M\t*\t0\t0\t*\t*\n"
            "readB\t0\tchr1\t250\t255\t50M\t*\t0\t0\t*\t*\n"
        )
        with open(self.sam_file, "w") as f:
            f.write(content)

    def read_sam_file(self):
        with open(self.sam_file, "r") as f:
            return f.readlines()

    def test_process_paired_sam_flags(self):
        process_paired_sam_flags(self.sam_file)
        result = self.read_sam_file()

        # Both mates of readA are unmapped, so each is also flagged
        # mate-unmapped (0x8): 1+4+8+64 = 77 and 1+4+8+128 = 141.
        # readB is mapped on both mates: 1+64 = 65 and 1+128 = 129.
        expected = [
            "@HD\tVN:1.0\tSO:unsorted\n",
            "@SQ\tSN:chr1\tLN:248956422\n",
            "readA\t77\tchr1\t100\t255\t50M\t*\t0\t0\t*\t*\n",
            "readA\t141\tchr1\t150\t255\t50M\t*\t0\t0\t*\t*\n",
            "readB\t65\tchr1\t200\t255\t50M\t*\t0\t0\t*\t*\n",
            "readB\t129\tchr1\t250\t255\t50M\t*\t0\t0\t*\t*\n",
        ]

        self.assertEqual(result, expected)


class TestProcessPairedSamFlagsEdgeCases(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.temp_sam_dir = tempfile.mkdtemp()
        self.sam_file = os.path.join(self.temp_sam_dir, "test.sam")

    def tearDown(self):
        shutil.rmtree(self.temp_sam_dir)
        super().tearDown()

    def _write_sam(self, records):
        header = "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"
        body = "".join(
            "\t".join(
                [name, str(flag)]
                + ["chr1", "100", "255", "50M", "*", "0", "0", "*", "*"]
            )
            + "\n"
            for name, flag in records
        )
        with open(self.sam_file, "w") as file:
            file.write(header + body)

    def _read_records(self):
        with open(self.sam_file) as file:
            return [
                line.rstrip("\n").split("\t")
                for line in file
                if not line.startswith("@")
            ]

    def test_orphaned_mate_is_dropped_without_error(self):
        # Filtering can remove one mate of a pair; the survivor has nothing to
        # be paired with, and must not take the next read's mate slot
        self._write_sam([("readA", 0), ("readB", 0), ("readB", 0)])
        process_paired_sam_flags(self.sam_file)
        records = self._read_records()

        self.assertEqual([record[0] for record in records], ["readB", "readB"])

    def test_odd_record_count_does_not_raise(self):
        # A lone trailing record used to be paired with an empty line
        self._write_sam([("readA", 0)])
        process_paired_sam_flags(self.sam_file)

        self.assertEqual(self._read_records(), [])

    def test_secondary_and_supplementary_records_are_dropped(self):
        # These are extra alignments of a read its primary record already
        # stands for, so keeping them would break the pairing
        self._write_sam([("readA", 0), ("readA", 256), ("readA", 2048), ("readA", 0)])
        process_paired_sam_flags(self.sam_file)
        records = self._read_records()

        self.assertEqual(len(records), 2)
        self.assertEqual([record[1] for record in records], ["65", "129"])

    def test_existing_mate_bits_are_not_doubled_up(self):
        # Minimap2 already sets first/second in pair; OR-ing the other bit on
        # top would mark a record as both, which samtools fastq discards
        self._write_sam([("readA", 97), ("readA", 145)])
        process_paired_sam_flags(self.sam_file)
        records = self._read_records()

        # Exactly one of the two mate bits must be set on each record
        for record in records:
            flag = int(record[1])
            self.assertTrue(bool(flag & 0x40) != bool(flag & 0x80))

    def test_returns_the_number_of_records_actually_written(self):
        # readB is an orphan: its partner was filtered out earlier, so it is
        # dropped here and must not be counted as retained
        self._write_sam([("readA", 0), ("readA", 0), ("readB", 0)])
        written = process_paired_sam_flags(self.sam_file)

        self.assertEqual(written, 2)
        self.assertEqual(len(self._read_records()), 2)

    def test_record_carrying_both_pair_bits_does_not_displace_its_partner(self):
        # Flag 192 sets first- and second-in-pair at once, so it satisfies both
        # selections; the pair must not collapse onto that one record
        self._write_sam([("readA", 192), ("readA", 0)])
        process_paired_sam_flags(self.sam_file)
        records = self._read_records()

        self.assertEqual(len(records), 2)
        first, second = records
        self.assertTrue(int(first[1]) & 0x40)
        self.assertTrue(int(second[1]) & 0x80)
        self.assertIsNot(first, second)

    def test_mates_are_ordered_by_their_pair_bits(self):
        # The second-in-pair record appears first in the file here; the output
        # must still put first-in-pair first
        self._write_sam([("readA", 145), ("readA", 97)])
        process_paired_sam_flags(self.sam_file)
        records = self._read_records()

        self.assertTrue(int(records[0][1]) & 0x40)
        self.assertTrue(int(records[1][1]) & 0x80)


class TestProcessSamFileSecondaryRecords(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        self.temp_sam_dir = tempfile.mkdtemp()
        self.sam_file = os.path.join(self.temp_sam_dir, "test.sam")

    def tearDown(self):
        shutil.rmtree(self.temp_sam_dir)
        super().tearDown()

    def test_unmapped_output_excludes_secondary_records(self):
        # A secondary alignment of a mapped read must not be emitted as an
        # extra copy of that read
        with open(self.sam_file, "w") as file:
            file.write("@HD\tVN:1.0\n")
            file.write("readA\t0\tchr1\t100\t60\t50M\t*\t0\t0\t*\t*\n")
            file.write("readA\t256\tchr1\t200\t0\t50M\t*\t0\t0\t*\t*\n")
            file.write("readB\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")

        process_sam_file(self.sam_file, "unmapped", None)

        with open(self.sam_file) as file:
            names = [line.split("\t")[0] for line in file if not line.startswith("@")]

        self.assertEqual(names, ["readB"])


if __name__ == "__main__":
    unittest.main()
