# LILAP
Low-input Low-cost amplification-free library-production method for PacBio Long-read sequencing

All the codes for LILAP are available in the 'code.sh' file. The installation of the software should be conducted with the correct versions specified in the accompanying paper.
The Snakemake file integrates automatic workflows for genome assembly, polishing, evaluation, SV calling, etc., with test data for ISO1-1 available in the snakemake/data directory. If errors occur, run commands individually, as they may stem from version differences of Snakemake or the programs used.


Links for the single-fly genome (and Wolbachia) assemblies:

ISO1-1: https://drive.google.com/file/d/1-BD99hcTxMuByG2PVFB0xjcttYHY1gKu/view?usp=drive_link

ISO1-2: https://drive.google.com/file/d/1Fi6L-yJmRprQ_virCgi18QCLBRF7HUZJ/view?usp=drive_link

wISO1-1: https://drive.google.com/file/d/1f2nLqUaaNm7orrgbCTvWwSh3s43uYKLr/view?usp=drive_link

wISO1-2: https://drive.google.com/file/d/1IjXDbNlbdbVfbI8CDNQcFazp7PgJY5Yb/view?usp=drive_link

All PacBio HiFi sequencing data generated in this study are available at the NCBI Sequence Read Archive database (https://www.ncbi.nlm.nih.gov/sra) under the accession numbers PRJNA983717.


## Usage
* Git clone all the codes, and run the following codes on the Linux server. The necessary software for this research is detailed in the "requirements.txt" file.
* Snakemake were prepared for genome assembly, evaluation, automatic polishing, structural variation detecting, etc.



## data for testing
Test the assembly code using both ISO1-1 and ISO1-2 ccs reads with the provided instructions.

For SVs and SNPs, utilize both the single-fly and Wolbachia genome assemblies. Examine the anticipated structural variation results in the sv.C01.ccs.hifiasm.txt file. It's important to note that some raw results necessitate additional manual verification, as elaborated in the manuscript. The complete pipeline, including the manual checking process, requires a few days for the analysis of the single-fly genome assembly.


## Outline
- [LILAP](#LILAP)
  * [code.sh](#codesh)
    + [ccs reads length distribution](#ccs-reads-length-distribution)
    + [Tn5 insertion site bias](#tn5-insertion-site-bias)
    + [ccs concensus](#ccs-concensus)
    + [relative depth](#relative-depth)
    + [GC bias](#gc-bias)
    + [assembly](#assembly)
      - [downsampling](#downsampling)
      - [hifiasm](#hifiasm)
      - [quast](#quast)
      - [BUSCO v5.4.2](#busco-v542)
      - [merqury calculate qv](#merqury-calculate-qv)
        * [calculate best kmer size](#calculate-best-kmer-size)
        * [build k-mer dbs with meryl](#build-k-mer-dbs-with-meryl)
        * [calculate qv by merqury](#calculate-qv-by-merqury)
    + [polish](#polish)
    + [Y chromosome identification](#y-chromosome-identification)
    + [SV calling](#sv-calling)
    + [non-B DNA detection](#non-b-dna-detection)
    + [SNP calling](#snp-calling)
      - [snpEff](#snpeff)
      - [get genome background bed](#get-genome-background-bed)
      - [verification and vaf filter of SNPs](#verification-and-vaf-filter-of-snps)
      - [MNP detecting](#mnp-detecting)
    + [Wolbachia & other bacteria](#wolbachia---other-bacteria)
      - [blobtools](#blobtools)





### ccs read length distribution
```sh
cat C01.subreads.cut19.fasta | awk '{if($0 ~ /^>/){rname=$0} else{all[rname]=all[rname]+length($0)}}END{for(rname in all){print all[rname]}}' > C01.subreads.cut19.fasta.length

cat C01.ccs.cut19.fasta | awk '{if($0 ~ /^>/){rname=$0} else{all[rname]=all[rname]+length($0)}}END{for(rname in all){print all[rname]}}' > C01.ccs.cut19.fasta.length
```


### Tn5 insertion site bias
```sh
perl Tn5_insertion_site_bias.pl
```


### ccs consensus
```sh
python3 ccs_circle.count.final.py C01.ccs.cut19.fasta C01.subreads.cut19.fasta C01

minimap2 -ax map-hifi --MD --eqx --secondary=no -t 10 C01.asm.p_ctg.fa C01.ccs.cut19.fasta | samtools view -F4 -F0x900 > C01.ccs.vs.asm.eqx.F4F0x900.sam

python3 rmdup.py C01

python3 ccs_identity.final.py C01.ccs.vs.asm.eqx.F4F0x900.rmdup.sam 
```


### relative depth
```sh
minimap2 -ax map-hifi --MD --secondary=no dm6.fa C01.ccs.cut19bp.fastq | samtools view -bS -F0x904 > C01.ccs.cut19bp.fq.dm6.F904.bam

samtools sort C01.ccs.cut19bp.fq.dm6.F904.bam -o C01.ccs.cut19bp.fq.dm6.F904.s.bam

samtools index C01.ccs.cut19bp.fq.dm6.F904.s.bam

samtools depth -a C01.ccs.cut19bp.fq.dm6.F904.s.bam > C01.ccs.cut19bp.fq.dm6.F904.s.depth

rm C01.ccs.cut19bp.fq.dm6.F904.bam


perl relative_depth_1.pl

perl relative_depth_2.pl

perl relative_depth_3.pl
```


### GC bias
```sh
perl GC_bias_1.pl C01.ccs.cut19bp.fq.dm6.F904.s.depth C01.ccs ccs

perl GC_bias_2.pl C01.ccs.500.dgc C01.ccs.500.dgc1 ccs
```


### assembly
#### downsampling
```sh
python3 downsampling.py C01.ccs.cut19.fasta 15
```


#### hifiasm
```sh
hifiasm -t40 -o C01.ccs.cut19bp.asm -f0 C01.ccs.cut19.fasta

awk '/^S/{print ">"$2;print $3}' C01.ccs.cut19bp.asm.p_ctg.gfa > C01.ccs.cut19bp.asm.p_ctg.fa

awk '/^S/{print ">"$2;print $3}' C01.ccs.cut19bp.asm.a_ctg.gfa > C01.ccs.cut19bp.asm.a_ctg.fa
```

#### quast
```sh
mkdir C01.ccs.cut19bp.asm

quast.py --large -t 40 -r dm6.fa -o C01.ccs.cut19bp.asm C01.ccs.cut19bp.asm.p_ctg.fa
```


#### BUSCO v5.4.2
```sh
mkdir C01.ccs.cut19bp.asm.pri.BUSCO542

busco -i C01.ccs.cut19bp.asm.p_ctg.fa -l diptera_odb10 -o C01.ccs.cut19bp.asm.pri.BUSCO542 -m genome -f
```


#### merqury calculate qv
##### calculate best kmer size
```sh
sh $MERQURY/best_k.sh 144000000
```


##### build k-mer dbs with meryl 
```sh
meryl count k=18 output C01.ccs.cut19bp.k18.meryl C01.ccs.cut19bp.fasta

meryl union-sum output C01.ccs.cut19bp.k18.merge.meryl C01.ccs.cut19bp.k18.meryl
```


##### calculate qv by merqury
```sh
merqury.sh C01.ccs.cut19bp.k18.merge.meryl C01.ccs.cut19bp.asm.p_ctg.fa C01.ccs.cut19bp.hifiasm
```


### polishing
```sh
meryl count k=15 C01.ccs.cut19bp.asm.p_ctg.fa output C01_lDB

meryl print greater-than distinct=0.9998 C01_lDB > repetitive_k15_C01_lDB.txt

winnowmap --MD -W repetitive_k15_C01_lDB.txt -ax map-pb -t20 C01.ccs.cut19bp.asm.p_ctg.fa C01.ccs.cut19.fasta > C01.ccs.cut19bp.pasm.sam

samtools sort -@20 -m2G -T C01.ccs.cut19bp.pasm.tmp -O bam -o C01.ccs.cut19bp.pasm.sort.bam C01.ccs.cut19bp.pasm.sam

rm C01.ccs.cut19bp.pasm.sam


falconc bam-filter-clipped -F=0x104 -t --output-count-fn=C01.ccs.cut19bp.pasm.sort.bam.filtered_aln_count.txt --output-fn=C01.ccs.cut19bp.pasm.sort.falconcF104.sam --input-fn=C01.ccs.cut19bp.pasm.sort.bam


/softwares/racon_liftover/build/bin/racon -t 20 \
C01.ccs.cut19.fasta \
C01.ccs.cut19bp.pasm.sort.falconcF104.sam \
C01.asm.p_ctg.fa \
-L C01.asm.p_ctg.racon.fa > C01.asm.p_ctg.racon.fa


jellyfish count -m 19 -t 10 -s 10000000000 -C C01.ccs.cut19.fasta -o C01.ccs.cut19.fasta.jf

jellyfish histo -t 10 C01.ccs.cut19.fasta.jf > C01.ccs.cut19.fasta.histo

genomescope2 -i C01.ccs.cut19.fasta.histo -o genomescope2_output -k 19


merfin -polish \
-sequence C01.asm.p_ctg.fa \
-seqmers C01.asm.p_ctg.meryl \
-readmers C01.ccs.cut19.meryl \
-peak 30.3  \
-vcf C01.asm.p_ctg.racon.fa.vcf \
-output C01.asm.p_ctg.fa.merfin.out \
-threads 20
```


### Y chromosome identification
```sh
cat C01.0.ccs.fasta C01.1.ccs.fasta C01.4.ccs.fasta > fly_6.ccs.male.fasta

cat C01.2.ccs.fasta C01.3.ccs.fasta C01.5.ccs.fasta > fly_6.ccs.female.fasta


minimap2 -ax map-hifi --MD --secondary=no -t 20 C01.asm.p_ctg.fa fly_6.ccs.male.fasta | samtools view -bS -F0x904 > C01_6.ccs.male.C01.contigs.F904.bam 

samtools sort C01_6.ccs.male.C01.contigs.F904.bam -o C01_6.ccs.male.C01.contigs.F904.s.bam 

rm C01_6.ccs.male.C01.contigs.F904.bam

samtools depth -a C01_6.ccs.male.C01.contigs.F904.s.bam  > C01_6.ccs.male.C01.contigs.depth


minimap2 -ax map-hifi --MD --secondary=no -t 20 C01.asm.p_ctg.fa fly_6.ccs.female.fasta | samtools view -bS -F0x904 > C01_6.ccs.female.C01.contigs.F904.bam 

samtools sort C01_6.ccs.female.C01.contigs.F904.bam -o C01_6.ccs.female.C01.contigs.F904.s.bam 

rm C01_6.ccs.female.C01.contigs.F904.bam

samtools depth -a C01_6.ccs.female.C01.contigs.F904.s.bam  > C01_6.ccs.female.C01.contigs.depth


python3 select_XY.median.py C01 contigs
```


### SV calling
```sh
nucmer --threads 40 --maxmatch --prefix dm62C01 dm6.fa C01.ccs.cut19bp.asm.p_ctg.fa

lastz dm6.fa[multiple] C01.ccs.cut19bp.asm.p_ctg.fa[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > dm62C01_lastz.txt

svmu dm62C01.delta dm6.fa C01.ccs.cut19bp.asm.p_ctg.fa h dm62C01_lastz.txt C01.ccs.hifiasm 

perl SV_results_filter.pl
```


### non-B DNA motif detection
```sh
cat TE_insertion_all.all.bed | awk '{print $1"\t"$2-1001"\t"$3+1000}' > TE_insertion_all.all.2000.bed

bedtools getfasta -fi /rd/yacai/dm6.fa -bed TE_insertion_all.all.2000.bed -fo TE_insertion_all.all.2000.fasta

gfa -seq TE_insertion_all.all.2000.fasta -out TE_insertion_all.all.2000.fasta.out

cat TE_insertion_all.all.2000.fasta.out_*.tsv | grep -v "Source" > TE_insertion_all.all.2000.fasta.all.out 


perl allmotifSum1.pl TE_insertion_all.all.2000 > TE_insertion_all.all.2000.fasta.xls1 

python3 non-B_DNA_sum2.py
```


### SNP calling
#### snpEff
```sh
java -jar snpEff.jar BDGP6.32.105 -i vcf pav_c01.txt -csvStats pav_c01.txt.ann.csv -stats pav_c01.txt.ann.html > C01.pav.ann.vcf
```


#### get genome background bed
```sh
cat rmsk.txt | awk '{ss=$7+1; ee=$8; if({if(($6=="chr2L" || $6=="chr2R" || $6=="chr3L" || $6=="chr3R" || $6=="chr4" || $6=="chrX" || $6=="chrY") && ($12=="Satellite")){print $6"\t"ss"\t"ee"\t"$12}}' | bedtools sort | bedtools merge | awk '{print $0"\tSatellite"}' > dm6.rmsk.Satellite.s.merge.bed

cat rmsk.txt | awk '{ss=$7+1; ee=$8; if({if(($6=="chr2L" || $6=="chr2R" || $6=="chr3L" || $6=="chr3R" || $6=="chr4" || $6=="chrX" || $6=="chrY") && ($12=="tSimple_repeat")){print $6"\t"ss"\t"ee"\t"$12}}' | bedtools sort | bedtools merge | awk '{print $0"\tSimple_repeat"}' > dm6.rmsk.Simple_repeat.s.merge.bed

cat rmsk.txt | grep -v "#" | awk '{if(($6=="chr2L" || $6=="chr2R" || $6=="chr3L" || $6=="chr3R" || $6=="chr4" || $6=="chrX" || $6=="chrY") && ($12=="LTR" || $12=="LINE" || $12=="SINE" || $12=="DNA" || $12 == "RC")){print $6"\t"$7"\t"$8"\tTE"}}' | awk '{ss=$2-50; ee=$3+50; print $1"\t"ss"\t"ee}' | bedtools sort | bedtools merge | awk '{print $0"\tTE"}' > dm6.TE.s.merge.50bp.s.merge.bed

cat simple_repeat.tsv | grep -v "#" | awk '{ss=$3+1; ee=$4; print $2"\t"ss"\t"ee"\tTRF"}' | bedtools sort | bedtools merge | awk '{print $0"\tSimple_repeat"}' > dm6.simple_repeat.trf.s.merge.bed

cat dm6.rmsk.Satellite.s.merge.bed dm6.rmsk.Simple_repeat.s.merge.bed dm6.TE.s.merge.50bp.s.merge.bed dm6.simple_repeat.trf.s.merge.bed | awk '{print $1"\t"$2"\t"$3}' | bedtools sort | bedtools merge > dm6.TE.s.merge.50bp.s.merge.simple_repeat.trf.s.merge.rmsk.Satellite_Simple_repeat.merge.bed
```


#### verification and vaf filter of SNPs
```sh
minimap2 -ax map-hifi --MD --secondary=no dm6.fa C01.ccs.cut19bp.fastq | samtools view -bS -F0x904 > C01.ccs.cut19bp.fq.dm6.F904.bam

samtools sort C01.ccs.cut19bp.fq.dm6.F904.bam -o C01.ccs.cut19bp.fq.dm6.F904.s.bam

samtools index C01.ccs.cut19bp.fq.dm6.F904.s.bam

rm C01.ccs.cut19bp.fq.dm6.F904.bam


python3 filter_sa_bam1.py C01.ccs.cut19bp.fq.dm6.F904.s.bam ### out file: C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.bam

python3 filter_sa_bam.py C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.bam ### out file: C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.rm_inside_te_simple.bam


cat C01.snp.out | awk '{print $1"\t"$2}' > C01.snp.pos

samtools mpileup --output-MQ --min-MQ 1 --min-BQ 0 -f dm6.fa -l C01.snp.pos C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.rm_inside_te_simple.bam > C01_E01.C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.rm_inside_te_simple.bam.baseq0.mapq0.pileup

python3 validated_non_rmsk_satellite_simple_trf_snp.baseq0.mapq0.py C01.snp.out C01_E01.C01 C01_E01.C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.rm_inside_te_simple.bam.baseq0.mapq0.pileup

python3 out_overlap_bed.py C01_E01.C01.filter.norm.merfin.filter.afdp.ann.vcf.shared_sp.snp.gene.repeat.repeattype.trf.rm_inside_te_simple_repeat.mapq0.baseq0.out dm6.rmsk.Satellite_Simple_repeat.trf.sm.up1500.bed up1500_vntr

cat C01_E01.C01.filter.norm.merfin.filter.afdp.ann.vcf.shared_sp.snp.gene.repeat.repeattype.trf.rm_inside_te_simple_repeat.mapq0.baseq0.up1500_vntr.out | awk '{if($18=="PASS" && $19=="out_up1500_vntr"){truetag="PASS"} else{truetag="FILTER"}}{print $0"\t"truetag}' > C01_E01.C01.filter.norm.merfin.filter.afdp.ann.vcf.shared_sp.snp.gene.repeat.repeattype.trf.rm_inside_te_simple_repeat.mapq0.baseq0.up1500_vntr.out1

python3 get_assembly_snp_vaf.baseq0.mapq0.py C01_E01.C01.snp.gene.repeat.repeattype.trf.rm_inside_te_simple_repeat.mapq0.baseq0.up1500_vntr.TRUE.rm_shared_diff.out C01_E01.C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.rm_inside_te_simple.bam.baseq0.mapq0.pileup

cat C01_E01.C01.filter.norm.merfin.filter.afdp.ann.vcf.shared_sp.snp.gene.repeat.repeattype.trf.rm_inside_te_simple_repeat.mapq0.baseq0.up1500_vntr.TRUE.vaf.out | awk '{if($23>2){print $0}}' > C01_E01.C01.shared_sp.snp.gene.repeat.repeattype.trf.rm_inside_te_simple_repeat.mapq0.baseq0.up1500_vntr.TRUE.vaf.minreads3.out
```


#### MNP detecting
```sh
samtools mpileup --output-MQ  --min-MQ 1 --min-BQ 0 --output-QNAME -f dm6.fa -l C01_E01.pav.min3.C01_sp.pos C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.rm_inside_te_simple.bam > C01_E01.pav.min3.C01_sp.pos.qname.pileup

for j in 1-9 1-99 1-299 1-1000
do
    filename="C01_E01.pav.min3.C01_sp.gt.out"
    pileup_file="C01_E01.pav.min3.C01_sp.pos.qname.pileup"
    python3 same_reads_phasing.py $filename $pileup_file $j
done

python3 combine_mnp_files.same_reads.py C01_sp
```


### Wolbachia & other bacteria
```sh
minimap2 -t 10 --secondary=no -ax map-hifi GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19.fasta | samtools view -S -F4 | awk '{print ">"$1; print $10}' > C01.ccs.cut19bp.wol.fasta

minimap2 -t 10 --secondary=no -ax asm5 GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.p_ctg.fa | samtools view -S -F4 | awk '{print ">"$1; print $10}' > C01.ccs.cut19bp.asm.wol.fa


nucmer --threads 40 --maxmatch --prefix wMel2C01asm GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.wol.fa

lastz GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.wol.fa --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > wMel2C01asm_lastz.txt

svmu wMel2C01asm.delta GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.wol.fa h wMel2C01asm_lastz.txt C01.asm.wol

dnadiff -p wMel2C01asm -d wMel2C01asm.delta
```


#### blobtools
```sh
blastn -db /blobtoolkit/nt \
       -query C01.ccs.cut19.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out C01.ncbi.blastn.out

diamond blastx \
        --query C01.ccs.cut19.fasta \
        --db /blobtoolkit/uniprot \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 16 \
        > C01.diamond.blastx.out

blobtools create \
    --fasta C01.ccs.cut19.fasta \
    --meta C01create.yaml \
    C01

blobtools view --remote C01
```


## Citation

To cite this repo in publications, please use

> Unpublished


## Contact

For questions or support, please contact caiyingao21@ioz.ac.cn.
