# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd

from q2_minimap2.types._format import (
    Minimap2IndexDBFmt,
    PairwiseAlignmentMN2DirectoryFormat,
    PairwiseAlignmentMN2Format,
)

from .test_minimap2 import MinimapTestsBase


class TestPaf(MinimapTestsBase):
    def test_PairwiseAlignmentMN2_Format_validate(self):
        filepath = self.get_data_path("type/paf_file.paf")
        format = PairwiseAlignmentMN2Format(filepath, mode="r")
        format.validate()

    def test_PairwiseAlignmentMN2DirectoryFormat_to_df_transformer(self):
        transformer = self.get_transformer(
            PairwiseAlignmentMN2DirectoryFormat, pd.DataFrame
        )
        paf = PairwiseAlignmentMN2DirectoryFormat(self.get_data_path("type/"), "r")
        paf_df = transformer(paf)
        self.assertIsInstance(paf_df, pd.DataFrame)

    def test_df_to_PairwiseAlignmentMN2DirectoryFormat_transformer(self):
        filepath = self.get_data_path("type/paf_file.paf")
        transformer = self.get_transformer(
            pd.DataFrame, PairwiseAlignmentMN2DirectoryFormat
        )
        paf_df = pd.read_csv(filepath)
        paf_format = transformer(paf_df)
        self.assertIsInstance(paf_format, PairwiseAlignmentMN2DirectoryFormat)


class TestIndex(MinimapTestsBase):
    def test_Minimap2IndexDBFmt_Format_validate(self):
        filepath = self.get_data_path("type/index.mmi")
        format = Minimap2IndexDBFmt(filepath, mode="r")
        format.validate()
