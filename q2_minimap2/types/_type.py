# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import FeatureData
from qiime2.plugin import SemanticType

Minimap2IndexDB = SemanticType("Minimap2IndexDB")
PairwiseAlignmentMN2 = SemanticType(
    "PairwiseAlignmentMN2", variant_of=FeatureData.field["type"]
)
