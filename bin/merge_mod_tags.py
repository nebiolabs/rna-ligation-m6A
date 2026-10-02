#!/usr/bin/env python3

################################################################################
# Illumina and ONT sequencing of RNA ligation products
# 
# Copyright (C) 2026 New England Biolabs
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import pysam
import argparse


def restore_modification_tags(original_bam, aligned_bam, output_bam):
    """
    Restore MM and ML tags from original BAM to aligned BAM based on read names.
    
    Args:
        original_bam (str): Path to original BAM with MM/ML tags.
        aligned_bam (str): Path to newly aligned BAM.
        output_bam (str): Path to output BAM with restored tags.
    """
    # Open BAM files
    orig = pysam.AlignmentFile(original_bam, "rb", check_sq=False)
    aligned = pysam.AlignmentFile(aligned_bam, "rb")
    out = pysam.AlignmentFile(output_bam, "wb", header=aligned.header)

    # Build dictionary of MM/ML tags from original BAM
    tag_dict = {}
    for read in orig.fetch(until_eof=True):
        tags = dict(read.tags)
        mm_tag = tags.get("MM")
        ml_tag = tags.get("ML")
        if mm_tag or ml_tag:
            tag_dict[read.query_name] = {"MM": mm_tag, "ML": ml_tag}

    # Iterate through aligned BAM and restore tags
    restored_count = 0
    for read in aligned.fetch(until_eof=True):
        if read.query_name in tag_dict:
            tags_to_restore = tag_dict[read.query_name]
            for tag, value in tags_to_restore.items():
                if value is not None:
                    read.set_tag(tag, value)
            restored_count += 1
        out.write(read)

    # Close files
    orig.close()
    aligned.close()
    out.close()

    print(f"Restored tags for {restored_count} reads. Output written to {output_bam}")


def main():
    parser = argparse.ArgumentParser(
        description="Restore MM and ML modification tags from original BAM to aligned BAM"
    )
    parser.add_argument("original_bam", help="Path to original BAM with MM/ML tags")
    parser.add_argument("aligned_bam", help="Path to newly aligned BAM")
    parser.add_argument("output_bam", help="Path to output BAM with restored tags")
    
    args = parser.parse_args()
    
    restore_modification_tags(args.original_bam, args.aligned_bam, args.output_bam)


if __name__ == "__main__":
    main()
