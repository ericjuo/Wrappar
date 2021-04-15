
rule migec_run:
    input:
        unpack(get_fastq)
    output:
        "pipeline/{sample}/{time}/migec/Assemble/MIG{size}/assemble.cmd.txt"
    params:
        checkout   = "pipeline/{sample}/{time}/migec/Checkout",
        hist       = "pipeline/{sample}/{time}/migec/Histogram",
        assem      = "pipeline/{sample}/{time}/migec/Assemble/MIG{size}",
        migec      = config["migec"],
        barcodes   = config["barcodes"],
        mig_size   = config["mig_size"],
        assem_meta = config["assembly_meta"]
    threads: 4
    resources:
         time   = 60,
         mem_mb = 18000
    shell:
        '''
            set -e
            
            java -jar -Xmx16G {params.migec} \
               Checkout \
               -cute {params.barcodes} \
               {input} {params.checkout}

           java -jar -Xmx16G {params.migec} \
               Histogram \
               {params.checkout} {params.hist}

           java -jar -Xmx16G {params.migec} \
               AssembleBatch \
               --force-collision-filter \
               --force-overseq {params.mig_size} \
               --sample-metadata {params.assem_meta} \
               {params.checkout} {params.hist} {params.assem}
        '''
 
