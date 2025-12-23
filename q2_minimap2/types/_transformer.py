# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import pandas as pd

from ..plugin_setup import plugin
from ._format import PairwiseAlignmentMN2DirectoryFormat


@plugin.register_transformer
def _1(data: pd.DataFrame) -> PairwiseAlignmentMN2DirectoryFormat:
    # Create a new PairwiseAlignmentMN2DirectoryFormat object
    ff = PairwiseAlignmentMN2DirectoryFormat()
    PairwiseAlignmentMN2_file_path = os.path.join(str(ff), "output.paf")
    data.to_csv(PairwiseAlignmentMN2_file_path, sep="\t", index=False, header=False)
    return ff


@plugin.register_transformer
def _2(data: PairwiseAlignmentMN2DirectoryFormat) -> pd.DataFrame:
    file = pd.read_csv(os.path.join(str(data), "output.paf"), sep="\t", header=None)
    return file
