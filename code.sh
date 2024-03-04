#!/bin/bash

### ccs reads length distribution
cat C01.subreads.cut19.fasta | awk '{if($0 ~ /^>/){rname=$0} else{all[rname]=all[rname]+length($0)}}END{for(rname in all){print all[rname]}}' > C01.subreads.cut19.fasta.length
cat C01.ccs.cut19.fasta | awk '{if($0 ~ /^>/){rname=$0} else{all[rname]=all[rname]+length($0)}}END{for(rname in all){print all[rname]}}' > C01.ccs.cut19.fasta.length


### Tn5 insertion site bias
perl Perl-10.pl


### ccs concensus
python3 ccs_circle.count.final.py C01.ccs.cut19.fasta C01.subreads.cut19.fasta C01
minimap2 -ax map-hifi --MD --eqx --secondary=no -t 10 C01.asm.p_ctg.fa C01.ccs.cut19.fasta | samtools view -F4 -F0x900 > C01.ccs.vs.asm.eqx.F4F0x900.sam
python3 rmdup.py C01
python3 ccs_identity.final.py C01.ccs.vs.asm.eqx.F4F0x900.rmdup.sam 


### relative depth
minimap2 -ax map-hifi --MD --secondary=no dm6.fa C01.ccs.cut19bp.fastq | samtools view -bS -F0x904 > C01.ccs.cut19bp.fq.dm6.F904.bam
samtools sort C01.ccs.cut19bp.fq.dm6.F904.bam -o C01.ccs.cut19bp.fq.dm6.F904.s.bam
samtools index C01.ccs.cut19bp.fq.dm6.F904.s.bam
samtools depth -a C01.ccs.cut19bp.fq.dm6.F904.s.bam > C01.ccs.cut19bp.fq.dm6.F904.s.depth
rm C01.ccs.cut19bp.fq.dm6.F904.bam

perl Perl-2.pl
perl Perl-13.pl
perl Perl-12.pl


### GC bias
perl Perl-10_3.pl C01.ccs.cut19bp.fq.dm6.F904.s.depth C01.ccs ccs
perl Perl-10_4.pl C01.ccs.500.dgc C01.ccs.500.dgc1 ccs


### assembly
#### downsampling
python3 downsampling.py C01.ccs.cut19.fasta 15

#### hifiasm
hifiasm -t40 -o C01.ccs.cut19bp.asm -f0 C01.ccs.cut19.fasta
awk '/^S/{print ">"$2;print $3}' C01.ccs.cut19bp.asm.p_ctg.gfa > C01.ccs.cut19bp.asm.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' C01.ccs.cut19bp.asm.a_ctg.gfa > C01.ccs.cut19bp.asm.a_ctg.fa

#### quast
mkdir C01.ccs.cut19bp.asm
quast.py --large -t 40 -r dm6.fa -o C01.ccs.cut19bp.asm C01.ccs.cut19bp.asm.p_ctg.fa

#### BUSCO v5.4.2
mkdir C01.ccs.cut19bp.asm.pri.BUSCO542
busco -i C01.ccs.cut19bp.asm.p_ctg.fa -l diptera_odb10 -o C01.ccs.cut19bp.asm.pri.BUSCO542 -m genome -f

#### merqury calculate qv
##### calculate best kmer size
sh $MERQURY/best_k.sh 144000000

##### build k-mer dbs with meryl 
meryl count k=18 output C01.ccs.cut19bp.k18.meryl C01.ccs.cut19bp.fasta
meryl union-sum output C01.ccs.cut19bp.k18.merge.meryl C01.ccs.cut19bp.k18.meryl

##### calculate qv by merqury
merqury.sh C01.ccs.cut19bp.k18.merge.meryl C01.ccs.cut19bp.asm.p_ctg.fa C01.ccs.cut19bp.hifiasm


### polish
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
-L C01.asm.p_ctg.racon.fa \
> C01.asm.p_ctg.racon.fa

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


### Y chromosome identification
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



### SV calling
nucmer --threads 40 --maxmatch --prefix dm62C01 dm6.fa C01.ccs.cut19bp.asm.p_ctg.fa
lastz dm6.fa[multiple] C01.ccs.cut19bp.asm.p_ctg.fa[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > dm62C01_lastz.txt
svmu dm62C01.delta dm6.fa C01.ccs.cut19bp.asm.p_ctg.fa h dm62C01_lastz.txt C01.ccs.hifiasm 
perl Perl-14.pl


### non-B DNA detection
cat TE_insertion_all.all.bed | awk '{print $1"\t"$2-1001"\t"$3+1000}' > TE_insertion_all.all.2000.bed
bedtools getfasta -fi /rd/yacai/dm6.fa -bed TE_insertion_all.all.2000.bed -fo TE_insertion_all.all.2000.fasta
gfa -seq TE_insertion_all.all.2000.fasta -out TE_insertion_all.all.2000.fasta.out
cat TE_insertion_all.all.2000.fasta.out_*.tsv | grep -v "Source" > TE_insertion_all.all.2000.fasta.all.out 

perl allmotifSum1.pl TE_insertion_all.all.2000 > TE_insertion_all.all.2000.fasta.xls1 
python3 non-B_DNA_sum2.py

### SNP calling
#### pav

#### deepvariant


#### snpEff
java -jar snpEff.jar BDGP6.32.105 -i vcf pav_c01.txt -csvStats pav_c01.txt.ann.csv -stats pav_c01.txt.ann.html > C01.pav.ann.vcf



#### get genome background bed
cat rmsk.txt | awk '{ss=$7+1; ee=$8; if({if(($6=="chr2L" || $6=="chr2R" || $6=="chr3L" || $6=="chr3R" || $6=="chr4" || $6=="chrX" || $6=="chrY") && ($12=="Satellite")){print $6"\t"ss"\t"ee"\t"$12}}' | bedtools sort | bedtools merge | awk '{print $0"\tSatellite"}' > dm6.rmsk.Satellite.s.merge.bed
cat rmsk.txt | awk '{ss=$7+1; ee=$8; if({if(($6=="chr2L" || $6=="chr2R" || $6=="chr3L" || $6=="chr3R" || $6=="chr4" || $6=="chrX" || $6=="chrY") && ($12=="tSimple_repeat")){print $6"\t"ss"\t"ee"\t"$12}}' | bedtools sort | bedtools merge | awk '{print $0"\tSimple_repeat"}' > dm6.rmsk.Simple_repeat.s.merge.bed
cat rmsk.txt | grep -v "#" | awk '{if(($6=="chr2L" || $6=="chr2R" || $6=="chr3L" || $6=="chr3R" || $6=="chr4" || $6=="chrX" || $6=="chrY") && ($12=="LTR" || $12=="LINE" || $12=="SINE" || $12=="DNA" || $12 == "RC")){print $6"\t"$7"\t"$8"\tTE"}}' | awk '{ss=$2-50; ee=$3+50; print $1"\t"ss"\t"ee}' | bedtools sort | bedtools merge | awk '{print $0"\tTE"}' > dm6.TE.s.merge.50bp.s.merge.bed
cat simple_repeat.tsv | grep -v "#" | awk '{ss=$3+1; ee=$4; print $2"\t"ss"\t"ee"\tTRF"}' | bedtools sort | bedtools merge | awk '{print $0"\tSimple_repeat"}' > dm6.simple_repeat.trf.s.merge.bed
cat dm6.rmsk.Satellite.s.merge.bed dm6.rmsk.Simple_repeat.s.merge.bed dm6.TE.s.merge.50bp.s.merge.bed dm6.simple_repeat.trf.s.merge.bed | awk '{print $1"\t"$2"\t"$3}' | bedtools sort | bedtools merge > dm6.TE.s.merge.50bp.s.merge.simple_repeat.trf.s.merge.rmsk.Satellite_Simple_repeat.merge.bed

#### verification and vaf filter of SNPs
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


#### MNP detecting
samtools mpileup --output-MQ  --min-MQ 1 --min-BQ 0 --output-QNAME -f dm6.fa -l C01_E01.pav.min3.C01_sp.pos C01.ccs.cut19bp.fq.dm6.F904.s.Fsa.rm_inside_te_simple.bam > C01_E01.pav.min3.C01_sp.pos.qname.pileup
for j in 1-9 1-99 1-299 1-1000
do
    filename="C01_E01.pav.min3.C01_sp.gt.out"
    pileup_file="C01_E01.pav.min3.C01_sp.pos.qname.pileup"
    python3 same_reads_phasing.py $filename $pileup_file $j
done
python3 combine_mnp_files.same_reads.py C01_sp


### Wolbachia & other bacteria
minimap2 -t 10 --secondary=no -ax map-hifi GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19.fasta | samtools view -S -F4 | awk '{print ">"$1; print $10}' > C01.ccs.cut19bp.wol.fasta
minimap2 -t 10 --secondary=no -ax asm5 GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.p_ctg.fa | samtools view -S -F4 | awk '{print ">"$1; print $10}' > C01.ccs.cut19bp.asm.wol.fa

nucmer --threads 40 --maxmatch --prefix wMel2C01asm GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.wol.fa
lastz GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.wol.fa --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > wMel2C01asm_lastz.txt
svmu wMel2C01asm.delta GCF_016584425.1_ASM1658442v1_genomic.fna C01.ccs.cut19bp.asm.wol.fa h wMel2C01asm_lastz.txt C01.asm.wol
dnadiff -p wMel2C01asm -d wMel2C01asm.delta

#### blobtools
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

./blobtools2/blobtools create \
    --fasta C01.ccs.cut19.fasta \
    --meta C01create.yaml \
    C01

./blobtools2/blobtools view --remote C01

