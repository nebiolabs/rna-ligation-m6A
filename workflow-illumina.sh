#!/bin/bash

################################################################################
# Two-way ligation (reaction 1)
################################################################################

apptainer exec apptainer.sif \
    minimap2 -a -x sr \
    references/way2.fasta data/two_way_ligation1.fastq.gz > results/two_way_ligation1.sam

apptainer exec apptainer.sif \
    samtools view -Sb results/two_way_ligation1.sam \
    | apptainer exec apptainer.sif \
        samtools sort > results/two_way_ligation1.bam

apptainer exec apptainer.sif \
    bin/context.py \
        --ctx-start 6 \
        --ctx-end 56 \
        --cut-site 31 \
        results/two_way_ligation1.bam references/way2.fasta > results/two_way_ligation1.tsv


################################################################################
# Two-way ligation (reaction 2)
################################################################################

apptainer exec apptainer.sif \
    minimap2 -a -x sr \
    references/way2.fasta data/two_way_ligation2.fastq.gz > results/two_way_ligation2.sam

apptainer exec apptainer.sif \
    samtools view -Sb results/two_way_ligation2.sam \
    | apptainer exec apptainer.sif \
        samtools sort > results/two_way_ligation2.bam

apptainer exec apptainer.sif \
    bin/context.py \
        --ctx-start 6 \
        --ctx-end 56 \
        --cut-site 31 \
        results/two_way_ligation2.bam references/way2.fasta > results/two_way_ligation2.tsv


################################################################################
# Three-way ligation
################################################################################

apptainer exec apptainer.sif \
    minimap2 -a -x sr \
    references/way3.fasta data/three_way_ligation.fastq.gz > results/three_way_ligation.sam

apptainer exec apptainer.sif \
    samtools view -Sb results/three_way_ligation.sam \
    | apptainer exec apptainer.sif \
        samtools sort > results/three_way_ligation.bam

apptainer exec apptainer.sif \
    bin/context.py \
        --ctx-start 65 \
        --ctx-end 116 \
        --cut-site 75 --cut-site 106 \
        results/three_way_ligation.bam references/way3.fasta > results/three_way_ligation.tsv


################################################################################
# IVT control ligation
################################################################################

apptainer exec apptainer.sif \
    minimap2 -a -x sr \
    references/wayI.fasta data/ivt_control.fastq.gz > results/ivt_control.sam

apptainer exec apptainer.sif \
    samtools view -Sb results/ivt_control.sam \
    | apptainer exec apptainer.sif \
        samtools sort > results/ivt_control.bam

apptainer exec apptainer.sif \
    bin/context.py \
        --ctx-start 50 \
        --ctx-end 100 \
        --cut-site 75 \
        results/ivt_control.bam references/wayI.fasta > results/ivt_control.tsv
