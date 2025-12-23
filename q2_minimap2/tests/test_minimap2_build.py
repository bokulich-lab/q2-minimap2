# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from unittest.mock import MagicMock, patch

from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2.build_index import build_index, create_idx_build_cmd

from .test_minimap2 import MinimapTestsBase


class TestMinimap2Build(MinimapTestsBase):
    def setUp(self):
        super().setUp()
        # Retrieve the file path for the test DNA sequences
        self.genome = DNAFASTAFormat(
            self.get_data_path("minimap2_build/dna-sequences.fasta"), mode="r"
        )

    @patch("subprocess.run")
    def test_build(self, mocked_run_command):
        build_index(self.genome)
        mocked_run_command.assert_called_once()

    @patch("subprocess.run")
    def test_build_hifi_preset(self, mocked_run_command):
        build_index(self.genome, preset="map-hifi")
        mocked_run_command.assert_called_once()

    @patch("subprocess.run")
    def test_build_map_pb_preset(self, mocked_run_command):
        build_index(self.genome, preset="map-pb")
        mocked_run_command.assert_called_once()

    @patch("subprocess.run")
    def test_build_sr_preset(self, mocked_run_command):
        build_index(self.genome, preset="sr")
        mocked_run_command.assert_called_once()

    @patch("subprocess.run")
    def test_build_index_subprocess_failure(self, mocked_run_command):
        # Configuring the mock to raise a CalledProcessError
        mocked_run_command.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd="minimap2 -d ..."
        )

        # Calling the method under test
        with self.assertRaisesRegex(
            Exception,
            "An error was encountered while running main.nf",
        ):
            build_index(self.genome)

    def test_create_idx_build_cmd(self):
        # Mock the database and sequences
        mock_database = MagicMock()
        mock_database.path = MagicMock()
        mock_database.path.__truediv__.return_value = "/mock/path/to/database/index.mmi"

        mock_sequences = "/mock/path/to/sequences.fasta"
        preset = "map-ont"

        expected_cmd = [
            "minimap2",
            "-x",
            preset,
            "-d",
            "/mock/path/to/database/index.mmi",
            mock_sequences,
        ]

        # Call the function
        cmd = create_idx_build_cmd(mock_database, mock_sequences, preset)

        # Assert the command matches the expected command
        self.assertEqual(cmd, expected_cmd)
