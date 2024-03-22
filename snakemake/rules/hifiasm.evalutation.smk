rule calculate_fasta_length_ccs:
    input:
        reads="data/{sample}.ccs.cut19bp.fasta"
    output:
        length="results/summary/{sample}.ccs.fasta.readlength"
    shell:
        """
        cat {input.reads} | awk '{{if($$0 ~ /^>/){{rname=$$0}} else{{all[rname]=all[rname]+length($$0)}}}}END{{for(rname in all){{print all[rname]}}}}' > {output.length}
        """

rule count_ccs_circles:
    input:
        reads="data/{sample}.ccs.cut19bp.fasta",
        subreads="data/{sample}.subreads.fasta.names"
    output:
        zmw_count="results/summary/{sample}.zmw.count.out",
        ccs_circle="results/summary/{sample}.ccs.circle.count.out"
    shell:
        """
        python3 scripts/ccs_circle.count.final.py {input.subreads} {input.reads} results/summary/{wildcards.sample}
        """

rule align_ccs_to_assembly:
    input:
        assembly="results/asm/{sample}.ccs.asm.p_ctg.fa",
        reads="data/{sample}.ccs.cut19bp.fasta"
    output:
        sam="results/summary/{sample}.ccs.vs.asm.eqx.F905.sam"
    threads: 10
    shell:
        """
        minimap2 -ax map-hifi --MD --eqx --secondary=no -t {threads} {input.assembly} {input.reads} | \
        samtools view -F0x904 > {output.sam}
        """

rule calculate_ccs_identity:
    input:
        sam="results/summary/{sample}.ccs.vs.asm.eqx.F905.sam",
        zmw_count="results/summary/{sample}.zmw.count.out"
    output:
        "results/summary/{sample}.zmw.count.out.asm.mapq_more_0.out1"
    shell:
        """
        python3 scripts/ccs_identity.final.py {input.sam} {input.zmw_count} asm
        """

rule hifiasm_assembly:
    input:
        reads="data/{sample}.ccs.cut19bp.fasta"
    output:
        p_ctg_gfa="results/asm/{sample}.ccs.asm.p_ctg.gfa"
    params:
        asm_prefix="results/asm/{sample}.ccs.asm"
    shell:
        """
        hifiasm -t 20 -o {params.asm_prefix} -f0 {input.reads}
        """

rule convert_p_ctg_to_fasta:
    input:
        gfa = "results/asm/{sample}.ccs.asm.p_ctg.gfa"
    output:
        fa = "results/asm/{sample}.ccs.asm.p_ctg.fa"
    shell:
        """
        awk '/^S/{{print ">"$2;print $3}}' {input.gfa} > {output.fa}
        """

rule quast_evaluation:
    input:
        assembly="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        quast_report="results/asm/quast/{sample}.ccs.quast_LG/report.tsv"
    params:
        reference="data/dm6.fa",  # Adjust the path to your reference genome as necessary
        output_dir="results/asm/quast/{sample}.ccs.quast_LG"
    threads: 40
    shell:
        """
        mkdir -p {params.output_dir}
        quast.py --large -t {threads} -r {params.reference} -o {params.output_dir} {input.assembly}
        """

# /rd/caiya/softwares/busco-5.4.2/bin/busco
rule busco_evaluation:
    input:
        assembly="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        "results/asm/{sample}.ccs.BUSCO/short_summary.specific.diptera_odb10.{sample}.ccs.BUSCO.json"        
    params:
        lineage="diptera_odb10",  # Adjust as necessary for your specific organism group
        output_dir="results/asm/{sample}.ccs.BUSCO"  # Specify the output directory for BUSCO
    threads: 40
    shell:
        """
        mkdir -p {params.output_dir}
        busco -i {input.assembly} -l {params.lineage} -o {params.output_dir} -m genome -f
        """


rule meryl_count:
    input:
        reads="data/{sample}.ccs.cut19bp.fasta"
    output:
        "results/asm/{sample}.ccs.merqury/{sample}.ccs.k18.meryl/0x000000.merylData"
    params:
        kmer_db="results/asm/{sample}.ccs.merqury/{sample}.ccs.k18.meryl"
    shell:
        """
        meryl count k=18 output {params.kmer_db} {input.reads}
        """

rule merqury_evaluation:
    input:
        input_merylData="results/asm/{sample}.ccs.merqury/{sample}.ccs.k18.meryl/0x000000.merylData", 
        p_fa = "results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        qv="results/asm/{sample}.ccs.merqury/{sample}.ccs.hifiasm.p.qv"
    params:
        infa="../{sample}.ccs.asm.p_ctg.fa", 
        kmer_db="results/asm/{sample}.ccs.merqury/{sample}.ccs.k18.meryl",
        output_name="{sample}.ccs.hifiasm.p"  # Adjust as necessary
    shell:
        """
        cd results/asm/{wildcards.sample}.ccs.merqury
        merqury.sh {params.kmer_db} {params.infa} {params.output_name}
        """


