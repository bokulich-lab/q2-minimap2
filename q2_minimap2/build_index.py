# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

from q2_types.feature_data import DNAFASTAFormat

from q2_minimap2._utils import run_command
from q2_minimap2.types._format import Minimap2IndexDBDirFmt


def create_idx_build_cmd(database, sequences, preset):
    build_cmd = [
        "minimap2",
        "-x",
        preset,
        "-d",  # Option for generating index
        str(database.path / "index.mmi"),  # Output index path
        str(sequences),  # Input sequences file path
    ]

    return build_cmd


def build_index(
    reference: DNAFASTAFormat,
    preset: str = "map-ont",  # Preset to apply multiple options
) -> Minimap2IndexDBDirFmt:
    # Initialize a directory format object to store the Minimap2 index
    database = Minimap2IndexDBDirFmt()

    # Construct the command to build the Minimap2 index file
    build_cmd = create_idx_build_cmd(database, reference, preset)

    try:
        # Execute the command to create the Minimap2 index database
        run_command(build_cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while using Minimap2, "
            f"(return code {e.returncode}), please re-run with --verbose "
            "to see the Minimap2 stdout and stderr."
        )

    return database
