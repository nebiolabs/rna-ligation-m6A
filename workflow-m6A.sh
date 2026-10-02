#!/bin/bash

# cd ont_data

# ont1="/data/potapov/projects/251125.ont_m6A_calling_LizhiL/jobs-m6A_DRACH/00_demux_with_seqtagger"
# ont2="/data/potapov/projects/260106.ont_m6A_calling_LizhiL/jobs-m6A_DRACH/00_demux_with_seqtagger"

# ln -sf $ont1/bam_demux.bc_1.bam m6A_005_percent.bam
# ln -sf $ont1/bam_demux.bc_2.bam m6A_030_percent.bam
# ln -sf $ont1/bam_demux.bc_3.bam m6A_080_percent.bam

# ln -sf $ont2/bam_demux.bc_1.bam m6A_000_percent.bam
# ln -sf $ont2/bam_demux.bc_2.bam m6A_050_percent.bam
# ln -sf $ont2/bam_demux.bc_3.bam m6A_100_percent.bam

ont_dir="/data/potapov/projects/260309.rna_ligation_manuscript_data_upload_Lizhi/github_repo/ont_data"
reference="/data/potapov/projects/260309.rna_ligation_manuscript_data_upload_Lizhi/github_repo/references/HERC3_m6A_dna.fasta"

for ix in 005 030 050 080 100
do
    echo ""
    echo "Processing barcode ${ix} ..."

    # Convert unaligned BAM reads to FASTQ
    apptainer exec apptainer.sif \
        samtools fastq --threads 20 $ont_dir/m6A_${ix}_percent.bam > results/m6A_${ix}_percent.fastq

    # Align FASTQ reads to the reference using minimap2
    apptainer exec apptainer.sif \
        minimap2 -t 20 -a -x map-ont $reference results/m6A_${ix}_percent.fastq \
        | apptainer exec apptainer.sif \
            samtools view -Sb > results/m6A_${ix}_percent.bam

    # Transfer modifications tags (MM and ML) to aligned BAM files
    apptainer exec --bind $PWD:$PWD apptainer.sif \
        bin/merge_mod_tags.py \
            $ont_dir/m6A_${ix}_percent.bam results/m6A_${ix}_percent.bam results/m6A_${ix}_percent_mods.bam

    # Sort and index the final BAM file
    apptainer exec apptainer.sif \
        samtools sort results/m6A_${ix}_percent_mods.bam > results/m6A_${ix}_percent_mods.sorted.bam

    apptainer exec apptainer.sif \
        samtools index results/m6A_${ix}_percent_mods.sorted.bam

    # detect m6A modifications and generate stats and plots
    apptainer exec --bind $PWD:$PWD apptainer.sif \
        bin/process_mods_m6A.py \
            --model m6A_DRACH \
            --min-read-length 900 \
            --output results/m6A_${ix}_percent.m6A_DRACH.csv \
            --save-mods results/m6A_${ix}_percent.m6A_DRACH.mods.csv \
            $reference results/m6A_${ix}_percent_mods.sorted.bam
done
