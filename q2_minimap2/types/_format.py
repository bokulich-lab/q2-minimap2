# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin import ValidationError, model


class Minimap2IndexDBFmt(model.BinaryFileFormat):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _validate_(self, level):
        # It's not clear if there is any way to tell if a Minimap2 index is
        # correct or not.
        # Minimap2 does have an inspect method
        pass


Minimap2IndexDBDirFmt = model.SingleFileDirectoryFormat(
    "Minimap2IndexDBDirFmt", "index.mmi", Minimap2IndexDBFmt
)


class PairwiseAlignmentMN2Format(model.TextFileFormat):
    def _validate(self, n_records=None):
        with open(str(self), "r") as file:
            line_number = 0
            for line in file:
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
        self._validate()


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
