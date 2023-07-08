
rule mixcr_run:
    input:
        "pipeline/{sample}/{time}/mitools/MIG{size}/{chain}_merged.t1.cf.fastq"
    output:
        "pipeline/{sample}/{time}/mixcr/MIG{size}/{imgt_release}/{chain}/clones_contigs_{chain}.tsv",
    params:
        base    = "pipeline/{sample}/{time}/mixcr/MIG{size}/{imgt_release}/{chain}",
        mixcr   = config["mixcr"],
        lib_ver = config["imgt_release"]
    threads: 4
    resources:
         time   = 20,
         mem_mb = 18000
    shell:
        '''
            set -e
            
            {params.mixcr} align \
                --library {params.lib_ver} \
                -f \
                -s dog \
                --verbose \
                -p kAligner2 \
                -OmaxHits=1 \
                -OvParameters.geneFeatureToAlign=VTranscript \
                {input} {params.base}/align_{wildcards.chain}.vdjca \
                --report {params.base}/log_align_{wildcards.chain}.txt

            {params.mixcr} assemble \
                -f \
                {params.base}/align_{wildcards.chain}.vdjca \
                {params.base}/output_{wildcards.chain}.clna \
                --write-alignments \
                --report {params.base}/log_assemble_{wildcards.chain}.txt


            {params.mixcr} assembleContigs \
                -f \
                --report {params.base}/log_assemblecontigs_{wildcards.chain}.txt \
                -OsubCloningRegion=VDJRegion \
                {params.base}/output_{wildcards.chain}.clna \
                {params.base}/output_contigs_{wildcards.chain}.clns 

            {params.mixcr} exportClones \
                --chains {wildcards.chain} \
                --preset fullImputed \
                {params.base}/output_contigs_{wildcards.chain}.clns \
                {output[0]}

        '''

