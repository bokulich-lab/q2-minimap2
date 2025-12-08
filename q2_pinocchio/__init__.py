# ----------------------------------------------------------------------------
# Copyright (c) 2024, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .build_index import build_index
from .classify_consensus import _find_consensus_annotation, classify_consensus_minimap2
from .extract_reads import extract_reads
from .filter_reads import filter_reads
from .minimap2_search import minimap2_search

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"

__all__ = [
    "filter_reads",
    "build_index",
    "minimap2_search",
    "classify_consensus_minimap2",
    "_find_consensus_annotation",
    "extract_reads",
]
