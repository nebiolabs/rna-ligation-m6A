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

"""
Process modifications in BAM file for DRACH motifs.
"""

import argparse
import sys
import pysam
from collections import defaultdict
from tqdm import tqdm


def get_iupac_bases(code):
    """Return bases matching IUPAC ambiguity code."""
    iupac = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],  # puRine
        'Y': ['C', 'T'],  # pYrimidine
        'S': ['G', 'C'],  # Strong
        'W': ['A', 'T'],  # Weak
        'K': ['G', 'T'],  # Keto
        'M': ['A', 'C'],  # aMino
        'B': ['C', 'G', 'T'],  # not A
        'D': ['A', 'G', 'T'],  # not C
        'H': ['A', 'C', 'T'],  # not G
        'V': ['A', 'C', 'G'],  # not T
        'N': ['A', 'C', 'G', 'T']
    }
    return iupac.get(code.upper(), [])


def is_drach_motif(seq, pos):
    """
    Check if position is the center A of a DRACH motif.
    DRACH: D=[AGT], R=[AG], A, C, H=[ACT]
    """
    if pos < 2 or pos >= len(seq) - 2:
        return False
    
    if seq[pos] != 'A':
        return False
    
    # Check D at pos-2 (not C)
    if seq[pos - 2] not in get_iupac_bases('D'):
        return False
    
    # Check R at pos-1 (purine: A or G)
    if seq[pos - 1] not in get_iupac_bases('R'):
        return False
    
    # Check C at pos+1
    if seq[pos + 1] != 'C':
        return False
    
    # Check H at pos+2 (not G)
    if seq[pos + 2] not in get_iupac_bases('H'):
        return False
    
    return True


def find_drach_motifs(fasta_path):
    """Find all DRACH motifs in reference sequences."""
    drach_positions = defaultdict(set)
    
    with pysam.FastaFile(fasta_path) as fasta:
        for chrom in fasta.references:
            seq = fasta.fetch(chrom).upper()
            
            for i in range(len(seq)):
                if is_drach_motif(seq, i):
                    drach_positions[chrom].add(i)
                    # Report motif with 5-base context
                    motif_seq = seq[i-2:i+3]
                    print(f"{chrom}\t{i}\t{motif_seq}")
            
            print(f"Found {len(drach_positions[chrom])} DRACH motifs in {chrom}", 
                  file=sys.stderr)
    
    return drach_positions


def find_all_a_positions(fasta_path):
    """Find all A positions in reference sequences."""
    a_positions = defaultdict(set)
    
    with pysam.FastaFile(fasta_path) as fasta:
        for chrom in fasta.references:
            seq = fasta.fetch(chrom).upper()
            
            for i in range(len(seq)):
                if seq[i] == 'A':
                    a_positions[chrom].add(i)
            
            print(f"Found {len(a_positions[chrom])} A positions in {chrom}", 
                  file=sys.stderr)
    
    return a_positions


def process_bam_modifications(bam_path, fasta_path, target_positions, threshold, num_reads=None, min_read_length=None, save_mods_path=None):
    """Process BAM file and extract A+a? modifications in target positions."""
    # Initialize statistics dictionaries
    position_coverage = defaultdict(int)  # Count of aligned bases
    position_modified = defaultdict(int)  # Count of modifications above threshold
    position_composition = defaultdict(lambda: defaultdict(int))  # Base composition at positions
    
    # Open file for saving individual modifications if requested
    mods_file = None
    if save_mods_path:
        mods_file = open(save_mods_path, 'w')
        mods_file.write("Chromosome,Position,ModScore\n")
    
    print("# Chromosome\tPosition\tRefBase\tQueryPos\tQueryBase\tModType\tScore", file=sys.stdout)
    
    reads_processed = 0
    
    with pysam.AlignmentFile(bam_path, "rb") as bam, \
         pysam.FastaFile(fasta_path) as fasta:
        
        # Get total number of reads for progress bar
        total_reads = num_reads if num_reads is not None else bam.count(read_callback='all')
        
        # Create progress bar
        pbar = tqdm(total=total_reads, desc="Processing reads", unit="reads", file=sys.stderr)
        
        for read in bam.fetch():
            # Check if we've reached the read limit
            if num_reads is not None and reads_processed >= num_reads:
                break
            
            if read.is_unmapped:
                pbar.update(1)
                continue
            
            if read.is_secondary or read.is_supplementary:
                pbar.update(1)
                continue

            if read.is_reverse:
                pbar.update(1)
                continue
            
            if read.mapping_quality < 60:
                pbar.update(1)
                continue
            
            # check that read maps to the reference within first 10 bases
            if read.reference_start > 10:
                pbar.update(1)
                continue

            # check that read maps to the reference within last 10 bases
            if read.reference_end is not None and read.reference_end < (fasta.get_reference_length(read.reference_name) - 10):
                pbar.update(1)
                continue
            
            if min_read_length is not None and read.query_length < min_read_length:
                pbar.update(1)
                continue

            reads_processed += 1

            # Get modified bases (forward strand)
            if read.modified_bases_forward:

                # Get reference positions for this read
                ref_positions = read.get_reference_positions(full_length=True)
                chrom = read.reference_name

                # # Count coverage for all target positions covered by this read
                # for query_pos, ref_pos in enumerate(ref_positions):
                #     position_composition[(chrom, ref_pos)][read.query_sequence[query_pos]] += 1
                
                # # Count coverage for all target positions covered by this read
                # for query_pos, ref_pos in enumerate(ref_positions):
                #     if ref_pos is not None and chrom in target_positions:
                #         if ref_pos in target_positions[chrom]:
                #             position_key = (chrom, ref_pos)
                #             position_coverage[position_key] += 1
            
                # Look for A modifications (m6A is typically 'A+a?' or 'A+a')
                for mod_key, modifications in read.modified_bases_forward.items():
                    # mod_key is tuple like ('A', 0, 'a') where:
                    # - first element is the canonical base
                    # - third element is the modification code
                    canonical_base, _, mod_code = mod_key

                    if canonical_base == 'A' and mod_code == 'a':

                        # Count coverage for all target positions covered by this read
                        for query_pos, ref_pos in enumerate(ref_positions):
                            position_composition[(chrom, ref_pos)][read.query_sequence[query_pos]] += 1

                        # modifications is list of tuples (query_pos, mod_score)
                        for query_pos, mod_score in modifications:
                            # Convert query position to reference position
                            ref_pos = ref_positions[query_pos]

                            if ref_pos is not None:
                                # Only process if position is in target positions
                                if chrom in target_positions and ref_pos in target_positions[chrom]:
                                    position_key = (chrom, ref_pos)
                                    
                                    # Save individual modification if requested
                                    if mods_file:
                                        mods_file.write(f"{chrom},{ref_pos},{mod_score}\n")
                                    
                                    # Count modifications above threshold
                                    if mod_score >= threshold:
                                        position_modified[position_key] += 1

                                    position_coverage[position_key] += 1
            
            pbar.update(1)
        
        pbar.close()
    
    # Close modifications file if opened
    if mods_file:
        mods_file.close()
    
    print(f"Processed {reads_processed} reads", file=sys.stderr)
    
    return position_coverage, position_modified, position_composition


def write_statistics(output_path, target_positions, position_coverage, position_modified, position_composition):
    """Write statistics to CSV file."""
    with open(output_path, 'w') as f:
        f.write("Chromosome,Position,Coverage,Modified,Fraction,A_count,C_count,G_count,T_count,N_count\n")
        
        # Iterate through all target positions
        for chrom in sorted(target_positions.keys()):
            for pos in sorted(target_positions[chrom]):
                position_key = (chrom, pos)
                coverage = position_coverage.get(position_key, 0)
                modified = position_modified.get(position_key, 0)
                fraction = (modified / coverage) if coverage > 0 else 0.0
                base_counts = position_composition.get(position_key, {})
                a_count = base_counts.get('A', 0)
                c_count = base_counts.get('C', 0)
                g_count = base_counts.get('G', 0)
                t_count = base_counts.get('T', 0)
                n_count = base_counts.get('N', 0)
                f.write(f"{chrom},{pos},{coverage},{modified},{fraction},{a_count},{c_count},{g_count},{t_count},{n_count}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Process modifications in BAM file for DRACH motifs'
    )
    parser.add_argument('fasta', help='Reference sequence in FASTA format')
    parser.add_argument('bam', help='Aligned BAM file')
    parser.add_argument('--output', '-o', required=True, 
                        help='Output CSV file for statistics')
    parser.add_argument('--model', '-m', choices=['m6A_DRACH', 'm6A'], default='m6A_DRACH',
                        help='Model to use: m6A_DRACH (A in DRACH context) or m6A (all A positions) (default: m6A_DRACH)')
    parser.add_argument('--threshold', '-t', type=int, default=200,
                        help='Modification score threshold (default: 200)')
    parser.add_argument('--num-reads', '-n', type=int, default=None,
                        help='Number of reads to process (default: all reads)')
    # minimal read length could be added if needed
    parser.add_argument('--min-read-length', type=int, default=None,
                        help='Minimum read length to process (default: no minimum)')
    parser.add_argument('--save-mods', type=str, default=None,
                        help='Output CSV file for individual modifications')
    
    args = parser.parse_args()
    
    # Load reference and find target positions based on model
    print(f"Loading reference from {args.fasta}...", file=sys.stderr)
    print(f"Using model: {args.model}", file=sys.stderr)
    
    if args.model == 'm6A_DRACH':
        print("# Chromosome\tPosition\tMotif", file=sys.stderr)
        target_positions = find_drach_motifs(args.fasta)
        total_positions = sum(len(positions) for positions in target_positions.values())
        print(f"Total DRACH motifs found: {total_positions}", file=sys.stderr)
    else:  # m6A
        target_positions = find_all_a_positions(args.fasta)
        total_positions = sum(len(positions) for positions in target_positions.values())
        print(f"Total A positions found: {total_positions}", file=sys.stderr)
    
    # Process BAM file for modifications
    print(f"\nProcessing BAM file: {args.bam}", file=sys.stderr)
    print(f"Using modification threshold: {args.threshold}", file=sys.stderr)
    if args.num_reads is not None:
        print(f"Processing up to {args.num_reads} reads", file=sys.stderr)
    else:
        print(f"Processing all reads", file=sys.stderr)
    
    if args.save_mods:
        print(f"Saving individual modifications to: {args.save_mods}", file=sys.stderr)
    
    position_coverage, position_modified, position_composition = process_bam_modifications(
        args.bam, args.fasta, target_positions, args.threshold, args.num_reads, args.min_read_length, args.save_mods
    )
    
    # Write statistics
    print(f"\nWriting statistics to: {args.output}", file=sys.stderr)
    write_statistics(args.output, target_positions, position_coverage, position_modified, position_composition)
    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
