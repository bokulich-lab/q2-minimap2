# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase


class MinimapTestsBase(TestPluginBase):
    package = "q2_minimap2.tests"

    def setUp(self):
        super().setUp()
