rule calculate_fasta_length_ccs:
    input:
        reads="data/{sample}.ccs.fasta"
    output:
        length="results/summary/{sample}.ccs.fasta.readlength"
    shell:
        """
        cat {input.reads} | awk '{{if($$0 ~ /^>/){{rname=$$0}} else{{all[rname]=all[rname]+length($$0)}}}}END{{for(rname in all){{print all[rname]}}}}' > {output.length}
        """

rule count_ccs_circles:
    input:
        reads="data/{sample}.ccs.fasta",
        subreads="data/{sample}.subreads.fasta.names"
    output:
        zmw_count="results/summary/{sample}.zmw.count.out",
        ccs_circle="results/summary/{sample}.ccs.circle.count.out"
    shell:
        """
        python3 scripts/ccs_circle.count.final.py {input.subreads} {input.reads} {wildcards.sample}
        """

rule align_ccs_to_assembly:
    input:
        assembly="results/asm/{sample}.ccs.asm.p_ctg.fa",
        reads="data/{sample}.ccs.fasta"
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
        reads="data/{sample}.ccs.fasta"
    output:
        asm_prefix="results/asm/{sample}.ccs.asm", 
        p_ctg_gfa="results/asm/{sample}.ccs.asm.p_ctg.gfa"
    threads: 40
    shell:
        """
        hifiasm -t {threads} -o {output.asm_prefix} -f0 {input.reads}
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
        quast_report="results/asmquast/{sample}.ccs_quast_LG/report.tsv"
    params:
        reference="data/dm6.fa",  # Adjust the path to your reference genome as necessary
        output_dir="results/asmquast/{sample}.ccs_quast_LG"
    threads: 40
    shell:
        """
        mkdir -p {params.output_dir}
        quast.py --large -t {threads} -r {params.reference} -o {params.output_dir} {input.assembly}
        """

rule busco_evaluation:
    input:
        assembly="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        directory("results/asm/{sample}.ccs.BUSCO")
    params:
        lineage="diptera_odb10",  # Adjust as necessary for your specific organism group
        output_dir="results/asm/{sample}.ccs.BUSCO"  # Specify the output directory for BUSCO
    threads: 40
    shell:
        """
        mkdir -p {params.output_dir}
        busco -i {input.assembly} -l {params.lineage} -o {output} -m genome -f
        """

rule meryl_count:
    input:
        reads="data/{sample}.ccs.fasta"
    output:
        kmer_db="results/asm/{sample}.ccs.merqury/{sample}.ccs.k18.meryl"
    shell:
        """
        meryl count k=18 output {output.kmer_db} {input.reads}
        """


rule merqury_evaluation:
    input:
        kmer_db="results/asm/{sample}.ccs.merqury/{sample}.ccs.k18.meryl",
        p_fa = "results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        output_name="results/asm/{sample}.ccs.merqury/{sample}.ccs.hifiasm.p"  # Adjust as necessary
    shell:
        """
        merqury.sh {input.kmer_db} {input.p_fa} {output.output_name}
        """



