rule trim_input_fasta:
    input:
        fasta="data/{sample}.ccs.fasta"
    output:
        "data/{sample}.ccs.cut19bp.fasta"
    shell:
        "cat {input.fasta} | python3 scripts/cut_adapter_19bp.py > {output}"

