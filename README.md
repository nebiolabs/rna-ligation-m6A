# RNA ligation

## Illumina data analysis

1. To analyze RNA ligation fidelity, Illumina sequencing reads are mapped to the reference sequences corresponding to expected ligation products. FASTA sequences of expected ligartion products are proivded in `references` directory:

    ```
    references/
    ├── way2.fasta
    ├── way3.fasta
    └── wayI.fasta
    ```

2. For each expected junction, context sequences at both sides of the junction are extracted and tabulated. Included `workflow-illumina.sh` script includes commands used to generata resulting tables in `results/` directory:

    ```
    results/
    ├── ivt_control.tsv
    ├── three_way_ligation.tsv
    ├── two_way_ligation1.tsv
    └── two_way_ligation2.tsv
    ```

## ONT data analysis

1. ONT modification detection is done with [`dorado`](https://github.com/nanoporetech/dorado/) basecaller (version 1.3.1) using `m6A_DRACH` and `inosine_m6A` models:

```
# m6A_DRACH
dorado basecaller rna004_130bps_hac@v5.1.0 data/pod5 \
  --recursive \
  --modified-bases m6A_DRACH \
  --modified-bases-threshold 0.0 \
  --emit-sam \
  --output-dir dorado_basecalled_reads

# inosine_m6A
dorado basecaller rna004_130bps_hac@v5.1.0 data/pod5 \
  --recursive \
  --modified-bases inosine_m6A \
  --modified-bases-threshold 0.0 \
  --emit-sam \
  --output-dir dorado_basecalled_reads
```

2. Basecalled ONT reads are demultiplexed using [SeqTagger](https://github.com/novoalab/SeqTagger/) software.

3. Resulting demultiplexed BAM files with modification information encoded in standard `MM` and `ML` tags are processed as follows:
  
    a. Convert unialigned ONT BAM files to FASTQ files

    b. Align FASTQ reads to the reference sequence (`references/HERC3_m6A_dna.fasta`)

    c. Transfer modifications tags (MM and ML) to aligned BAM files

    d. Detect m6A modifications and generate stats

    e. Included `workflow-m6A.sh` script includes commands used to generata resulting tables in `results/` directory.
