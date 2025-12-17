# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pandas as pd
import pandas.testing as pdt
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from qiime2 import Artifact

from q2_pinocchio.classify_consensus import _PAF_format_df_to_series_of_lists

from .test_pinocchio import PinocchioTestsBase


def series_is_subset(expected, observed):
    joined = pd.concat([expected, observed], axis=1, join="inner")
    compared = joined.apply(lambda x: x.iloc[0].startswith(x.iloc[1]), axis=1)
    return len(compared[compared]) >= len(compared[~compared])


class TestConsensusAssignment(PinocchioTestsBase):
    def setUp(self):
        super().setUp()

        taxonomy_tsv_path = self.get_data_path("consensus/taxonomy.tsv")
        reads_path = self.get_data_path("consensus/se-dna-sequences.fasta")

        self.taxonomy = Artifact.import_data(FeatureData[Taxonomy], taxonomy_tsv_path)
        self.reads = Artifact.import_data(FeatureData[Sequence], reads_path)
        self.paf = Artifact.load(self.get_data_path("consensus/search_results.qza"))

    def test_classify_consensus_minimap2(self):
        (paf, taxonomy) = self.plugin.pipelines["classify_consensus_minimap2"](
            query=self.reads,
            reference=self.reads,
            reference_taxonomy=self.taxonomy,
        )

        expected = self.taxonomy.view(pd.Series)
        self.assertTrue(series_is_subset(expected, taxonomy.view(pd.Series)))

    def test_PAF_format_df_to_series_of_lists(self):
        paf = self.paf.view(pd.DataFrame)
        taxonomy = self.taxonomy.view(pd.Series)

        obs = _PAF_format_df_to_series_of_lists(paf, taxonomy)
        exp = pd.Series(
            {
                "1111561": [
                    "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; "
                    "o__Legionellales; f__; g__; s__",
                    "k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; "
                    "o__; f__; g__; s__",
                ],
                "835097": [
                    "k__Bacteria; p__Chloroflexi; c__SAR202; o__; f__; g__; s__"
                ],
                "junk": ["Unassigned"],
            }
        )

        # Use assert to compare the sorted lists
        assert exp.tolist() == obs.tolist(), "Contents of the Series are not equal."

    # should fail when hit IDs are missing from reference taxonomy
    def test_PAF_format_df_to_series_of_lists_fail_on_missing_ids(
        self,
    ):
        paf = self.paf.view(pd.DataFrame)
        taxonomy = self.taxonomy.view(pd.Series)
        paf = paf.astype(object)
        paf.loc[3] = ["junk", "", "", "", "", "lost-id"] + [""] * 17
        with self.assertRaisesRegex(KeyError, "results do not match.*lost-id"):
            _PAF_format_df_to_series_of_lists(paf, taxonomy)

    def test_find_consensus_annotation(self):
        (consensus,) = self.plugin.methods["_find_consensus_annotation"](
            search_results=self.paf, reference_taxonomy=self.taxonomy
        )
        obs = consensus.view(pd.DataFrame)

        exp = pd.DataFrame(
            {
                "Taxon": {
                    "1111561": "k__Bacteria; p__Proteobacteria; "
                    "c__Gammaproteobacteria",
                    "835097": "k__Bacteria; p__Chloroflexi; c__SAR202; o__; f__; "
                    "g__; s__",
                    "junk": "Unassigned",
                },
                "Consensus": {"1111561": "1.0", "835097": "1.0", "junk": "1.0"},
            }
        )
        pdt.assert_frame_equal(exp, obs, check_names=False)
