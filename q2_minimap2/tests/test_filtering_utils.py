# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
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
from unittest.mock import MagicMock, call, patch

from q2_minimap2._filtering_utils import (
    build_filtered_out_dir,
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

from .test_minimap2 import Minimap2TestsBase


class TestBuildFilteredOutDir(Minimap2TestsBase):
    @patch("os.listdir")
    @patch("shutil.copy")
    def test_copy_files(self, mock_copy, mock_listdir):
        # Set up
        input_reads = type("test", (object,), {"path": "fake_input_path"})
        filtered_seqs = type("test", (object,), {"path": "fake_filtered_path"})
        mock_listdir.return_value = ["file1.txt", "file2.txt"]

        # Invoke
        build_filtered_out_dir(input_reads, filtered_seqs)

        # Assert
        expected_calls = [
            call(os.path.join("fake_filtered_path", "file1.txt"), "fake_input_path"),
            call(os.path.join("fake_filtered_path", "file2.txt"), "fake_input_path"),
        ]
        mock_copy.assert_has_calls(expected_calls, any_order=True)


class TestSetPenalties(Minimap2TestsBase):

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


class TestCalculateIdentity(Minimap2TestsBase):
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


class TestGetAlignmentLength(Minimap2TestsBase):
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


class TestProcessSamFile(Minimap2TestsBase):
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


class TestConvertToFasta(Minimap2TestsBase):
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


class TestConvertToFastq(Minimap2TestsBase):
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


class TestMakeMn2Cmd(Minimap2TestsBase):
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
        ]

        # Call the function and compare the result with the expected command
        actual_cmd = make_mn2_cmd(
            mapping_preset, index, n_threads, penalties, reads, None, samf_fp
        )
        self.assertEqual(actual_cmd, expected_cmd)


class TestConvertToFastqPaired(Minimap2TestsBase):
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


class TestMakeMn2PairedEndCmd(Minimap2TestsBase):
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

        result_cmd = make_mn2_cmd(
            mapping_preset, index, n_threads, penalties, reads1, reads2, samf_fp
        )
        self.assertEqual(result_cmd, expected_cmd)


class TestRunCmd(Minimap2TestsBase):
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


class TestCollateSamInplace(Minimap2TestsBase):
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


class TestProcessPairedSamFlags(Minimap2TestsBase):
    def setUp(self):
        # Create a temporary SAM file for testing
        self.temp_dir = tempfile.mkdtemp()
        self.sam_file = os.path.join(self.temp_dir, "test.sam")
        self.create_sam_file()

    def tearDown(self):
        # Remove the temporary directory and file
        shutil.rmtree(self.temp_dir)

    def create_sam_file(self):
        content = (
            "@HD\tVN:1.0\tSO:unsorted\n"
            "@SQ\tSN:chr1\tLN:248956422\n"
            "read1\t4\tchr1\t100\t255\t50M\t*\t0\t0\t*\t*\n"
            "read2\t4\tchr1\t150\t255\t50M\t*\t0\t0\t*\t*\n"
            "read3\t0\tchr1\t200\t255\t50M\t*\t0\t0\t*\t*\n"
            "read4\t0\tchr1\t250\t255\t50M\t*\t0\t0\t*\t*\n"
        )
        with open(self.sam_file, "w") as f:
            f.write(content)

    def read_sam_file(self):
        with open(self.sam_file, "r") as f:
            return f.readlines()

    def test_process_paired_sam_flags(self):
        process_paired_sam_flags(self.sam_file)
        result = self.read_sam_file()

        expected = [
            "@HD\tVN:1.0\tSO:unsorted\n",
            "@SQ\tSN:chr1\tLN:248956422\n",
            "read1\t69\tchr1\t100\t255\t50M\t*\t0\t0\t*\t*\n",
            "read2\t133\tchr1\t150\t255\t50M\t*\t0\t0\t*\t*\n",
            "read3\t65\tchr1\t200\t255\t50M\t*\t0\t0\t*\t*\n",
            "read4\t129\tchr1\t250\t255\t50M\t*\t0\t0\t*\t*\n",
        ]

        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
