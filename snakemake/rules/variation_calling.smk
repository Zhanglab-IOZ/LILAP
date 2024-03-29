# Rule for nucmer
rule nucmer:
    input:
        ref="data/dm6.fa",
        query="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        delta="results/SV/dm62{sample}.delta"
    threads: 40
    params:
        prefix="results/SV/dm62{sample}"
    shell:
        "nucmer --threads {threads} --maxmatch --prefix {params.prefix} {input.ref} {input.query}"

# Rule for lastz
rule lastz:
    input:
        ref="data/dm6.fa",
        query="results/asm/{sample}.ccs.asm.p_ctg.fa"
    output:
        "results/SV/dm62{sample}_lastz.txt"
    shell:
        "lastz {input.ref}[multiple] {input.query}[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 > {output}"

# Rule for svmu
rule svmu:
    input:
        delta="results/SV/dm62{sample}.delta",
        ref="data/dm6.fa",
        query="results/asm/{sample}.ccs.asm.p_ctg.fa", 
        lastz="results/SV/dm62{sample}_lastz.txt"
    output:
        sv="results/SV/sv.{sample}.ccs.hifiasm.p.txt"
    params:
        prefix="{sample}.ccs.hifiasm.p",
        delta="dm62{sample}.delta",
        lastz="dm62{sample}_lastz.txt",
        ref="../../data/dm6.fa",
        query="../asm/{sample}.ccs.asm.p_ctg.fa"
    shell:
        """
        cd results/SV
        svmu {params.delta} {params.ref} {params.query} h {params.lastz} {params.prefix}
        """

# Rule for running the Perl script for SV results filtering
rule sv_results_filter:
    input:
        sv="results/SV/sv.{sample}.ccs.hifiasm.p.txt"
    output:
        svfilter="results/SV/sv.{sample}.ccs.hifiasm.p.filter0"  # Adjust output file name as needed
    shell:
        "perl scripts/SV_results_filter.pl {input.sv} > {output.svfilter}"
