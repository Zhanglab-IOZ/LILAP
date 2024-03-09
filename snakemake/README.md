
# Genome Assembly, Polishing, Evaluation and SV Detection Workflow

This repository contains a Snakemake workflow designed for genome assembly, genome evaluation, and structural variation (SV) detecting. It leverages state-of-the-art tools to provide a streamlined process for handling genomic data, from ccs reads to assembled genomes and their evaluation.

![Image text](https://github.com/Zhanglab-IOZ/LILAP/blob/2e3f13b794069753036001c1b48f7bd7037a98dc/snakemake/snakemake_workflow_graph.png)



## Overview

The workflow is divided into several key stages:

1. **Genome Assembly** - Utilizes tools like hifiasm for assembling genome sequences from PacBio Hi-Fi sequencing data (example in the data: C01.ccs.fasta).
2. **Sequencing Evaluation** - Some in-house python and perl scripts to access the ccs reads length distribution, relative sequencing depth, ccs reads QV, etc. 
3. **Genome Evaluation** - Employs various metrics and tools (quast for reference-based assembly evaluation, BUSCO for single copy orthologs completeness and merqury for reference-free QV evaluation) to assess the quality of the assembled genomes.
4. **Automatic Polishing** - A state-of-the-art automated polishing tool form Mc Cartney, A.M. 2022 et al. were used in the snakemake, which aim to correct the assembly errors and improve the QV of the genome assembly.
5. **SV Detection** - Identifies structural variations from the assembled genomes against the reference genome (example prepared in the data/dm6.fa) using a widely used program named svmu (Chakraborty, M. 2018 et al.).


Each stage is encapsulated in a separate Snakemake rule file, ensuring modularity and ease of use.

## Dependencies
We strongly encourage you to intall the dependencies via mamba instead of conda:
To install mamba, please fellow:
```
   conda install -c conda-forge mamba
```

- Snakemake v8.0.0
- Hifiasm v0.12
- Python v3.11.6
- Perl v5.32.1
- minimap2 v2.24-r1122
- samtools v1.17
- quast v5.2.0
- BUSCO v5.4.2
- meryl v2.3
- merqury v1.3
- winnowmap v2.03
- pb-falconc v1.13.1
- racon liftover branch v1.6.0
- jellyfish v2.2.10
- GenomeScope v2.0
- merfin v1.1
- bcftools v1.14
- mummer v4.0.0
- lastz v1.04.22
- svmu


## Installation

1. Clone this repository:
   ```
   git clone https://github.com/Zhanglab-IOZ/LILAP.git
   ```
2. Install Snakemake and other dependencies. You can install part of the dependencies via conda:
   ```
   conda create -n your_env_name -c bioconda --file requirement.txt
   ```
   However, some program can only download from the software author's github repository like svmu. If you find any questions for installation, please follow the guideline of the software official website.


## Usage

To run the entire workflow:

```
snakemake --cores N
```

If the snakemake workflow breaks unexpectedly (possibly due to a configuration of a dependencies, or some potential program bug, which happens quite often ╮(╯▽╰)╭ ), you can just keep running from the break point:

```
snakemake --cores N --rerun-incomplete
```

## Input Data

We prepare a demo data of ISO1-1 downsize data to test the entire pipeline, while you can download them from Google Drive:

   * D.melanogaster release 6 reference genome: [Google Drive Link](https://drive.google.com/file/d/1auUP206WUfA-Dba0Td-Fbr_1kdoAbt1M/view?usp=sharing)
   
   * C01.ccs.fasta: [Google Drive Link](https://drive.google.com/file/d/1hxBG3qVU1YBEDHVhvTciWc8GoN-eOCjO/view?usp=sharing)

   * C01.subreads.fasta.names (used for ccs identity analysis): [Google Drive Link](https://drive.google.com/file/d/1J7NSVweBkCzcdeTVnGGOUM9cw_hkT2D_/view?usp=sharing)

Please put the demo data in to the data/ directory

## Output
All the output file will located in the "results" directory:

Here are parts of the major results and their positions:
results
├── asm
│   ├── C01.ccs.asm.a_ctg.gfa # hifiasm alternative assembly result (.gfa file)
│   ├── C01.ccs.asm.p_ctg.fa  # hifiasm primary assembly result (.fa file)
│   ├── C01.ccs.asm.p_ctg.gfa # hifiasm primary assembly result (.gfa file)
│   ├── C01.ccs.BUSCO
│   │   ├── short_summary.specific.diptera_odb10.ultralow.asm.a_ctg.C01.ccs.BUSCO.json # BUSCO evaluatiuon result file of primary assembly (.json file)
│   │   └── short_summary.specific.diptera_odb10.ultralow.asm.a_ctg.C01.ccs.BUSCO.txt  # BUSCO evaluatiuon result file of primary assembly (.txt file)
│   ├── C01.ccs.merqury
│   │   ├── C01.ccs.hifiasm.p.completeness.stats # kmer completeness of the primary assembly from merqury
│   │   └── C01.ccs.hifiasm.p.qv # QV of the primary assembly from merqury
│   └── quast
│       └── C01.ccs.quast_LG
│           ├── report.html   # quastLG evaluation report (.html file)
│           ├── report.tex # quastLG evaluation report (.tex file)
│           ├── report.tsv # quastLG evaluation report (.tsv file)
│           └── report.txt # quastLG evaluation report (.txt file)
├── polish
│   ├── C01.asm.p_ctg.racon.fa   # automatic polishing intermediate file after racon (.fa file)
│   ├── C01.ccs.asm.p_ctg.fa.merfin.out.polish.fa  # automatic polishing final fasta, after racon & merfin (.fa file)
│   ├── C01.ccs.polish.BUSCO
│   │   ├── short_summary.specific.diptera_odb10.ultralow.asm.a_ctg.C01.ccs.polish.BUSCO.json # BUSCO evaluatiuon result file of final polishing fasta (.json file)
│   │   └── short_summary.specific.diptera_odb10.ultralow.asm.a_ctg.C01.ccs.polish.BUSCO.txt  # BUSCO evaluatiuon result file of final polishing fasta (.txt file)
│   ├── C01.ccs.polish.merqury
│   │   ├── C01.ccs.hifiasm.p.polish.completeness.stats  # kmer completeness of the final polishing fasta from merqury
│   │   └── C01.ccs.hifiasm.p.polish.qv   # QV of the final polishing fasta from merqury
│   └── quast
│       └── C01.ccs.polish.quast_LG
│           ├── report.html   # quastLG evaluation report (.html file)
│           ├── report.tex # quastLG evaluation report of final polishing fasta (.tex file)
│           ├── report.tsv # quastLG evaluation report of final polishing fasta (.tsv file)
│           └── report.txt # quastLG evaluation report of final polishing fasta (.txt file)
├── summary
│   ├── C01.ccs.circle.count.out # subreads pass count summary (column 1: subreads pass, column 2: count)
│   ├── C01.ccs.fasta.readlength # ccs reads read length distribution (read length, bp)
│   ├── C01.zmw.count.out  # zmw id subreads pass count summary (column 1: zmw id, column 2: subreads pass count)
│   └── C01.zmw.count.out.asm.mapq_more_0.out1  # ccs reads generated in a zmw self alignment against to primary assembly error rate (column 1: zmw id, column 2: subreads pass, column 3: identical base number, column 4: total read length)
└── SV
    ├── cm.C01.ccs.hifiasm.p.txt # A tab delimited file with all the reference genomic regions that are present in higher copy numbers (>1) in the sample genome. Those with "trans" in their names mean either it is a transposable element or non-TE copies of a gene in different chromosomes.
    ├── cnv_all.C01.ccs.hifiasm.p.txt  # A tab delimited file with all the reference genomic regions that are present in higher copy numbers (>1) in the sample genome. Those with "trans" in their names mean either it is a transposable element or non-TE copies of a gene in different chromosomes.
    ├── dm62C01.delta   # whole genome alignment delta file (from mummer)
    ├── dm62C01_lastz.txt  # whole genome alignment lastz file (from lastz)
    ├── small.C01.ccs.hifiasm.p.txt # A tab delimited file containing SNPs and small indels that occur within syntenic blocks (or MUMs).
    ├── sv.C01.ccs.hifiasm.p.filter0   # SV list after primary filter from the raw SV list
    └── sv.C01.ccs.hifiasm.p.txt # A tab delimited file that summarizes structural mutations (indels, CNVs, inversions) in the sample genome with respect to the reference genome.


## Citation

To cite this repo in publications, please use

> Unpublished


## Contact

For questions or support, please contact caiyingao21@ioz.ac.cn.
