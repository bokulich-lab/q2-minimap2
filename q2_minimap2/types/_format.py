# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2.plugin import ValidationError, model

# Number of mandatory (non-tag) columns in a PAF record
PAF_MANDATORY_FIELDS = 12

# Columns holding identifiers or symbols that must stay verbatim strings
PAF_STRING_COLUMNS = (0, 4, 5)


# Read a PAF file into a DataFrame.
# PAF records carry a variable number of trailing tag columns (a no-hit row has
# 12-13, a secondary alignment 22, a primary alignment 23). pandas fixes the
# column count from the first line and raises if a later line is wider, so the
# widest row has to be found before parsing. A run with no alignments at all
# leaves the file empty, which pandas also refuses to parse.
def read_paf(paf_fp):
    max_fields = 0
    with open(paf_fp, "r") as file:
        for line in file:
            if line.strip():
                max_fields = max(max_fields, line.count("\t") + 1)

    if max_fields == 0:
        return pd.DataFrame(columns=range(PAF_MANDATORY_FIELDS))

    # Sequence identifiers are read verbatim: left to pandas, a reference named
    # "NA" or "None" would become a missing value and one named "007" would
    # turn into the number 7, neither of which can be matched back to the
    # reference taxonomy.
    return pd.read_csv(
        paf_fp,
        sep="\t",
        header=None,
        names=range(max_fields),
        dtype={column: str for column in PAF_STRING_COLUMNS},
        na_filter=False,
    )


# Leading bytes of a Minimap2 index, the same magic number Minimap2 itself
# checks for before deciding whether a file is an index
MM_IDX_MAGIC = b"MMI\x02"


class Minimap2IndexDBFmt(model.BinaryFileFormat):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _validate_(self, level):
        # Handed a file that is not an index, Minimap2 silently falls back to
        # reading it as FASTA and reports every query as unmapped, so the magic
        # number is worth checking here rather than letting that surface as an
        # empty result. Reading four bytes costs the same at either level.
        with self.open() as file:
            magic = file.read(len(MM_IDX_MAGIC))

        if magic != MM_IDX_MAGIC:
            raise ValidationError(
                "File does not appear to be a Minimap2 index: expected it to "
                f"start with {MM_IDX_MAGIC!r}, found {magic!r}."
            )


Minimap2IndexDBDirFmt = model.SingleFileDirectoryFormat(
    "Minimap2IndexDBDirFmt", "index.mmi", Minimap2IndexDBFmt
)


class PairwiseAlignmentMN2Format(model.TextFileFormat):
    # How many records to inspect per validation level. A PAF from a long-read
    # run can hold millions of alignments, and QIIME asks for the "min" level
    # on every view, so scanning the whole file each time is what the level is
    # meant to avoid.
    _record_count_map = {"min": 5, "max": None}

    def _validate(self, n_records=None):
        with open(str(self), "r") as file:
            line_number = 0
            for line in file:
                if n_records is not None and line_number >= n_records:
                    break

                line_number += 1
                fields = line.strip().split("\t")

                # Check for at least 12 fields
                if len(fields) < 12:
                    raise ValidationError(
                        f"Line {line_number}: Insufficient number of fields."
                    )

                # Validate specific fields for type, range, and specific conditions
                try:
                    # Fields 2, 3, 4, 7, 8, 9, 10, 11 should be integers >= 0
                    query_seq_length = int(fields[1])  # Query sequence length
                    query_start = int(fields[2])  # Query start
                    query_end = int(fields[3])  # Query end
                    target_seq_length = int(fields[6])  # Target sequence length
                    target_start = int(fields[7])  # Target start
                    target_end = int(fields[8])  # Target end
                    matching_bases = int(fields[9])  # Number of matching bases
                    total_bases = int(fields[10])  # Number of bases, including gaps

                    # Ensure values are non-negative
                    for value in [
                        query_seq_length,
                        query_start,
                        query_end,
                        target_seq_length,
                        target_start,
                        target_end,
                        matching_bases,
                        total_bases,
                    ]:
                        if value < 0:
                            raise ValueError("Value cannot be negative.")

                    # Ensure query start is less than or equal to query end,
                    # and similarly for target
                    if query_start > query_end:
                        raise ValueError("Query start greater than query end.")
                    if target_start > target_end:
                        raise ValueError("Target start greater than target end.")

                    # Mapping quality must be an integer between 0 and 255
                    mq = int(fields[11])
                    if mq < 0 or mq > 255:
                        raise ValueError("Mapping quality out of bounds.")

                except ValueError as e:
                    raise ValidationError(f"Line {line_number}: {e}")

                # Check strand field to be '+' or '-' or '*'
                if fields[4] not in ["+", "-", "*"]:
                    raise ValidationError(
                        f'Line {line_number}: Strand field (5th column) must be "+" , '
                        f'"-" or "*" but is {fields[4]}.'
                    )

    def _validate_(self, level):
        self._validate(self._record_count_map[level])


# A directory format for PAF files where each file ends with .paf and
# is named according to the sample it represents.
class PairwiseAlignmentMN2DirectoryFormat(model.DirectoryFormat):
    PairwiseAlignmentMN2_files = model.FileCollection(
        r".+\.paf$", format=PairwiseAlignmentMN2Format
    )

    @PairwiseAlignmentMN2_files.set_path_maker
    def PairwiseAlignmentMN2_path_maker(self, sample_id):
        """
        Constructs a path for a PAF file using the provided sample_id.
        """
        return f"{sample_id}.paf"
