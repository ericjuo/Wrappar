
rule mitools_merge:
    input:
        assem_log = "pipeline/{sample}/{time}/migec/Assemble/MIG{size}/assemble.cmd.txt"
    output:
        "pipeline/{sample}/{time}/mitools/MIG{size}/{chain, (?!IGH).+}_merged.t1.cf.fastq"
    threads: 4
    resources:
         time   = 20,
         mem_mb = 12000
    run:
        import os

        base = os.path.dirname(input.assem_log)
        r2 = os.path.join(base,wildcards.chain+"_R2.t1.cf.fastq")
        r1 = r2.replace("_R2.t1.cf.fastq","_R1.t1.cf.fastq")

        if "IGH" not in wildcards.chain:
            shell(f'''
                java -Xmx16G -jar {config["mitools"]} merge -ss -s 0.7 \
                    {r1} {r2} {{output}}
            ''')


rule heavy_chain_merge:
    input:
        igmu = "pipeline/{sample}/{time}/mitools/MIG{size}/IGM_merged.t1.cf.fastq",
        igg  = "pipeline/{sample}/{time}/mitools/MIG{size}/IGG_merged.t1.cf.fastq",
        iga  = "pipeline/{sample}/{time}/mitools/MIG{size}/IGA_merged.t1.cf.fastq",
        ige  = "pipeline/{sample}/{time}/mitools/MIG{size}/IGE_merged.t1.cf.fastq"
    output:
        "pipeline/{sample}/{time}/mitools/MIG{size}/IGH_merged.t1.cf.fastq"
    shell:
        '''
            cat {input.igmu} {input.igg} {input.iga} {input.ige} > {output}
        '''
