rule meryl_count_1DB_fa:
    input:
        p_ctg_fa="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        "results/polish/{sample}_lDB/0x000000.merylData"
    params:
        "results/polish/{sample}_lDB"
    shell:
        "meryl count k=15 {input.p_ctg_fa} output {params}"

rule meryl_count_1DB_fasta:
    input:
        "results/polish/{sample}_lDB"
    output:
        "results/polish/repetitive_k15_{sample}_lDB.txt"
    shell:
        "meryl print greater-than distinct=0.9998 {input} > {output}"

rule meryl_count_fa:
    input:
        p_ctg_fa="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        fa_meryl="results/polish/meryl/{sample}.asm.p_ctg.meryl"
    shell:
        "meryl count k=19 output {output.fa_meryl} {input.p_ctg_fa}"

rule meryl_count_fasta:
    input:
        fasta="data/{sample}.ccs.fasta"
    output:
        fasta_meryl="results/polish/meryl/{sample}.ccs.cut19.meryl"
    shell:
        "meryl count k=19 output {output.fasta_meryl} {input.fasta}"

rule winnowmap:
    input:
        fa="results/asm/{sample}.ccs.asm.p_ctg.fa",
        txt="results/polish/repetitive_k15_{sample}_lDB.txt",
        fasta="data/{sample}.ccs.fasta"
    output:
        "results/polish/{sample}.ccs.cut19bp.pasm.sam"
    shell:
        "winnowmap --MD -W {input.txt} -ax map-pb -t20 {input.fa} {input.fasta} > {output}"

rule samtools_sort:
    input:
        "results/polish/{sample}.ccs.cut19bp.pasm.sam"
    output:
        temp("results/polish/{sample}.ccs.cut19bp.pasm.sort.bam")
    shell:
        "samtools sort -@20 -m2G -T {sample}.ccs.cut19bp.pasm.tmp -O bam -o {output} {input}"

rule falconc:
    input:
        "results/polish/{sample}.ccs.cut19bp.pasm.sam"
    output:
        count_txt="results/polish/{sample}.ccs.cut19bp.pasm.sort.bam.filtered_aln_count.txt",
        sam="results/polish/{sample}.ccs.cut19bp.pasm.sort.falconcF104.sam"
    shell:
        "falconc bam-filter-clipped -F=0x104 -t --output-count-fn={output.count_txt} --output-fn={output.sam} --input-fn={input}"

rule racon:
    input:
        fasta="data/{sample}.ccs.fasta",
        sam="results/polish/{sample}.ccs.cut19bp.pasm.sort.falconcF104.sam",
        fa="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        out_fa="results/polish/{sample}.asm.p_ctg.racon.fa", 
        out_vcf="results/polish/{sample}.asm.p_ctg.racon.fa.vcf"
    shell:
        "/path/to/racon_liftover/build/bin/racon -t 20 {input.fasta} {input.sam} {input.fa} -L > {output.out_fa}"

rule jellyfish_count:
    input:
        "data/{sample}.ccs.fasta"
    output:
        "results/polish/{sample}.ccs.fasta.jf"
    shell:
        "jellyfish count -m 19 -t 10 -s 10000000000 -C {input} -o {output}"

rule jellyfish_histo:
    input:
        "results/polish/{sample}.ccs.fasta.jf"
    output:
        "results/polish/{sample}.ccs.fasta.histo"
    shell:
        "jellyfish histo -t 10 {input} > {output}"


# Rule to run genomescope2 and capture kcov
rule genomescope2_and_capture_kcov:
    input:
        histo_file="results/polish/{sample}.ccs.fasta.histo"
    output:
        kcov="results/polish/genomescope2_output/{sample}.kcov.txt",
    shell:
        """
        genomescope2 -i {input.histo_file} -o results/polish/genomescope2_output -k 19 > results/polish/genomescope2_output/log.txt
        awk '{{if(NR == 3){{split($4, a, ":"); print a[2]}}}}' results/polish/genomescope2_output/log.txt > {output.kcov}
        """

# Rule to run merfin using kcov from genomescope2
rule merfin_polish:
    input:
        fa="results/asm/{sample}.ccs.asm.p_ctg.fa", 
        seqmers="results/polish/meryl/{sample}.asm.p_ctg.meryl",
        readmers="results/polish/meryl/{sample}.ccs.cut19.meryl",
        kcov="results/polish/genomescope2_output/{sample}.kcov.txt", 
        vcf="results/polish/{sample}.asm.p_ctg.racon.fa.vcf"
    output:
        output_prefix="results/polish/{sample}.asm.p_ctg.fa.merfin.out", 
        output_vcf="results/polish/{sample}.asm.p_ctg.fa.merfin.out.polish.vcf"
    shell:
        """
        kcov=$(cat {input.kcov})
        merfin -polish \
        -sequence {input.fa} \
        -seqmers {input.seqmers} \
        -readmers {input.readmers} \
        -peak $kcov  \
        -vcf {input.vcf} \
        -output {output} \
        -threads 20
        """

# Rule for compressing the VCF file with bcftools view
rule compress_vcf:
    input:
        "results/polish/{sample}.asm.p_ctg.fa.merfin.out.polish.vcf"
    output:
        "results/polish/{sample}.asm.p_ctg.fa.merfin.out.polish.vcf.gz"
    shell:
        "bcftools view -Oz {input} > {output}"

# Rule for indexing the compressed VCF file with bcftools index
rule index_vcf:
    input:
        "results/polish/{sample}.asm.p_ctg.fa.merfin.out.polish.vcf.gz"
    output:
        "results/polish/{sample}.asm.p_ctg.fa.merfin.out.polish.vcf.gz.csi"
    shell:
        "bcftools index {input}"

# Rule for creating a consensus sequence with bcftools consensus
rule generate_consensus:
    input:
        vcf="results/polish/{sample}.asm.p_ctg.fa.merfin.out.polish.vcf.gz",
        fa="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        "results/polish/{sample}.ccs.asm.p_ctg.fa.merfin.out.polish.fa"
    shell:
        "bcftools consensus {input.vcf} -f {input.fa} -H 1 > {output}"



