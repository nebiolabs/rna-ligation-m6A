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

import argparse
import re
import pysam


# Corrected read classification (2026-10-02).
# Categorization priority: mismatch > indel > indel_at_junction
#   mismatch          - at least one mismatch anywhere (indels may be present)
#   indel             - no mismatches, at least one indel outside of the junction
#   indel_at_junction - no mismatches, all indels are adjacent to the junction
def classify_context(encoded_seq):
    if encoded_seq.find('X') != -1:
        return('mismatch')

    if encoded_seq.find('-') == -1 and encoded_seq.find('+') == -1:
        return('correct')

    # every run of consecutive indels must be adjacent to a cut site '|'
    for m in re.finditer(r'[-+]+', encoded_seq):
        s, e = m.span()
        at_junction = (s > 0 and encoded_seq[s - 1] == '|') or (e < len(encoded_seq) and encoded_seq[e] == '|')
        if not at_junction:
            return('indel')

    return('indel_at_junction')


parser = argparse.ArgumentParser()
parser.add_argument('bam')
parser.add_argument('reference')
parser.add_argument('--cut-site', type=int, action='append', default=[])
parser.add_argument('--ctx-start', type=int, default=0)
parser.add_argument('--ctx-end', type=int, default=150)
args = parser.parse_args()

reference = None

with pysam.FastxFile(args.reference) as fasta:
    for entry in fasta:
        reference = entry.sequence.upper()
        break

total = 0
counts = {}

with pysam.AlignmentFile(args.bam) as bamfile:
    for read in bamfile:
        if read.is_unmapped:
            continue

        if read.is_secondary:
            continue

        if read.is_supplementary:
            continue
        
        if read.reference_start > args.ctx_start:
            continue

        if read.reference_end < args.ctx_end:
            continue
        
        s = None
        e = None
        
        aligned_pairs = read.get_aligned_pairs()

        status = 0
        cut_site = 0

        aligned_seq = ''
        encoded_seq = ''

        for (qp,rp) in aligned_pairs:
            if rp == args.ctx_start:
                status = 1
                s = qp

            if rp == args.ctx_end:
                status = 0
                e = qp
            
            if status == 0:
                continue
            
            if (rp != None) and (rp in args.cut_site) and aligned_seq[-1] != '|':
                aligned_seq += '|'
                encoded_seq += '|'
                cut_site = 1

            if qp != None and rp != None:
                aligned_seq += read.query_sequence[qp].upper()

                if read.query_sequence[qp] == reference[rp]:
                    # correct base
                    encoded_seq += '.'
                else:
                    # substitution
                    encoded_seq += 'X'

            elif qp == None and rp != None:
                # deletion
                aligned_seq += '-'
                encoded_seq += '-'
            elif qp != None and rp == None:
                # insertion
                aligned_seq += read.query_sequence[qp].lower()
                encoded_seq += '+'
            else:
                print('error: qp == None and rp = None')
                exit()
            
        ctx = (aligned_seq,encoded_seq)

        if ctx not in counts:
            counts[ctx] = 0
        
        counts[ctx] += 1
        total += 1

print('context', 'encoding', 'type', 'count', 'fraction', sep='\t')

for ctx in sorted(counts, key=lambda x: counts[x], reverse=True):
    print(ctx[0], ctx[1], classify_context(ctx[1]), counts[ctx], counts[ctx] / total, sep='\t')
