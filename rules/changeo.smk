# RETHINK THIS RULE SO WE DONT HAVE TO DOWNLOAD THAT FASTA
rule make_blast_db:
    output:
       #igblast_fasta = "share/igblast/fasta/imgt_dog_ig_VDJ.fasta",
        gapped_fasta  = "share/imgt/ig_VDJ.fasta"
    params:
        imgt_dir = "share/imgt",
        fasta_dir = "share/igblast/fasta",
        imgt_url = config["imgt_url"],
       #igblast_dir = config["igblast_dir"]
    log:
       #"Share/IgBlast/make_blast_db.done"
        "logs/make_blast_db/make_blastdb.log"
    shell:
        '''
            set -e

            # log stdout and stderr
            mkdir -p logs/{rule}
            exec &> >(tee {log})

            mkdir -p {params.fasta_dir}

            for i in {{V,D,J}}; do
            
                wget {params.imgt_url}/IGH$i.fasta -P {params.imgt_dir}
                wget {params.imgt_url}/IGK$i.fasta -P {params.imgt_dir} || true
                wget {params.imgt_url}/IGL$i.fasta -P {params.imgt_dir} || true
            
                cat {params.imgt_dir}/IGH$i.fasta \
                    {params.imgt_dir}/IGK$i.fasta \
                    {params.imgt_dir}/IGL$i.fasta > {params.imgt_dir}/ig_$i.fasta || true
                
                
                edit_imgt_file.pl {params.imgt_dir}/ig_$i.fasta \
                        > {params.fasta_dir}/imgt_dog_ig_$i.fasta &&
                
                makeblastdb -parse_seqids \
                    -dbtype nucl \
                    -in {params.fasta_dir}/imgt_dog_ig_$i.fasta \
                    -out share/igblast/database/imgt_dog_ig_$i
            
            done 

            cat {params.imgt_dir}/ig_*.fasta > {output.gapped_fasta}
        '''

#rule make_gapped_blast_db:
#    output:
#        igblast_fasta = "share/igblast_gapped/fasta/imgt_dog_ig_VDJ.gapped.fasta",
#        gapped_fasta  = "share/imgt_gapped/ig_vdj.gapped.fasta"
#    params:
#        canis_rep = config["canis_rep"],
#        imgt_dir  = "share/imgt_gapped",
#        fasta_dir = "share/igblast_gapped/fasta",
#    log:
#        "logs/make_gapped_blast_db/make_gapped_blastdb.log"
#    run:
#        import os
#        import fileinput
#       
#        # get gapped fastas for each chain type (e.g. IGHV, IGHD, IGHJ, etc.)
#        fas = {}
#        for root, dirs, files in os.walk(params.canis_rep):
#            # only wanted gapped ig fastas - no T-cells
#            if ("fastaNGaps" in root) and ("TR" not in root):
#                for f in files:
#                    if f.endswith(".fasta"):
#                        fa = os.path.join(root,f)
#                        key = fa.split("/")[-3]
#                        if key not in fas:
#                            fas[key] = defaultdict(list)
#                            fas[key] = [fa]
#                        else:
#                            fas[key].append(fa)
#
#        # set and create imgt_gapped dir
#        base = params.imgt_dir
#        os.makedirs(base, exist_ok=True)
#
#        # concatenate fastas for each chain type
#        for k,v in fas.items():
#            if k not in ["IGHC","IGKC","IGLC"]:
#                out = os.path.join(base, f"{k}.gapped.fasta")
#                with open(out, "w") as fout, fileinput.input(v) as fin:
#                    for line in fin:
#                        fout.write(line)
#            
#        # new out name to edit fasta for makeblastdb
#       #db_name = f"imgt_dog_{k}"
#       #imgt_fasta = f"imgt_dog_{k}.gapped.fasta"
#       #edit_out = os.path.join(base, f"imgt_dog_{k}.gapped.fasta")
#        shell('''
#            set -e
#            
#            for i in {{V,D,J}}; do
#            
#                cat {params.imgt_dir}/IGH$i.gapped.fasta \
        #                    {params.imgt_dir}/IGK$i.gapped.fasta \
        #                    {params.imgt_dir}/IGL$i.gapped.fasta \
        #                    > {params.imgt_dir}/ig_$i.gapped.fasta || true
#            
#                edit_imgt_file.pl {params.imgt_dir}/ig_$i.gapped.fasta \
        #                        > {params.fasta_dir}/imgt_dog_ig_$i.gapped.fasta
#
#                makeblastdb -parse_seqids \
        #                    -dbtype nucl \
        #                    -in {params.fasta_dir}/imgt_dog_ig_$i.gapped.fasta \
        #                    -out share/igblast_gapped/database/imgt_dog_ig_gapped_$i
#            done
#        ''')
#
#        # concatenate all edited fastas for rule output and gapped fastas
#        # for MakeDb ref input below
#        shell('''
#            cat {params.fasta_dir}/* > {output.igblast_fasta}
#            cat {params.imgt_dir}/ig_*.gapped.fasta > {output.gapped_fasta}
#        ''')

rule convert_to_fasta:
    input:
        fastq = "pipeline/{sample}/{time}/mitools/MIG1/{chain}_merged.t1.cf.fastq"
    output:
        fasta = "analysis/{imgt_release}/{chain}/IgBlastn/mod_clones/{sample}_{time}_{chain}.fasta"
    params:
        tmp_dir = "analysis/{imgt_release}/{chain}/IgBlastn/mod_clones/tmp"
    log:
        "logs/{imgt_release}/convert_to_fasta/{sample}_{time}_{chain}.log"
    run:
        # due to the nature of preprocessing data, the mitools output contains
        # fastq seqs with some identical names from merging chains into igh.
        # first convert to fasta format
        tmp_fasta = os.path.join(params.tmp_dir,
                f"{wildcards.sample}_{wildcards.time}_{wildcards.chain}.tmp.fasta")

        shell(f'''
            set -e

            # log stdout and stderr
            mkdir -p logs/{rule}
            exec &> >(tee {log})
            
            mkdir -p {{params.tmp_dir}}
            seqtk seq -a {{input.fastq}} > {tmp_fasta}
        ''')
        
        # second need to add a number to part of the name
        with open(tmp_fasta,"r") as f, open(output.fasta,"w") as out:
            i = 1
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    line = line.split(" ")
                    line[0] = f"{line[0]}_{i}"
                    print(" ".join(line), file=out)
                    i += 1     
                else:
                    print(line, file=out)

rule run_igblastn:
    input:
       #"share/igblast/fasta/imgt_dog_ig_VDJ.fasta",
        gapped_fasta  = "share/imgt/ig_VDJ.fasta",
        fasta         = "analysis/{imgt_release}/{chain}/IgBlastn/mod_clones/{sample}_{time}_{chain}.fasta"
    output:
        blastn_fmt7 = "analysis/{imgt_release}/{chain}/IgBlastn/{sample}_{time}_{chain}.fmt7"
    params:
        igdata = config["igdata"]["dir"],
        aux    = config["igdata"]["aux"],
        db_dir = "share/igblast/database"
    log:
        "logs/{imgt_release}/run_igblastn/{sample}_{time}_{chain}.log"
    shell:
        '''
            set -e

            # log stdout and stderr
            mkdir -p logs/{rule}
            exec &> >(tee {log})
            
            export IGDATA={params.igdata}
           
           igblastn \
                -germline_db_V {params.db_dir}/imgt_dog_ig_V \
                -germline_db_D {params.db_dir}/imgt_dog_ig_D \
                -germline_db_J {params.db_dir}/imgt_dog_ig_J \
                -organism rhesus_monkey \
                -auxiliary_data {params.aux} \
                -domain_system imgt \
                -ig_seqtype Ig \
                -outfmt '7 std qseq sseq btop' \
                -query {input.fasta} \
                -num_alignments_V 1 \
                -num_alignments_D 1 \
                -num_alignments_J 1 \
                -out {output.blastn_fmt7}
        '''

rule changeo_makedb:
    input:
        gapped_fasta  = "share/imgt/ig_VDJ.fasta",
       #germline = "share/igblast/fasta/imgt_dog_ig_VDJ.fasta",
        mod_clones    = "analysis/{imgt_release}/{chain}/IgBlastn/mod_clones/{sample}_{time}_{chain}.fasta",
        blastn_fmt7   = "analysis/{imgt_release}/{chain}/IgBlastn/{sample}_{time}_{chain}.fmt7"
    output:
        "analysis/{imgt_release}/{chain}/MakeDb/{sample}_{time}_{chain}_db-pass.tsv"
    params:
        out_dir = "analysis/{imgt_release}/{chain}/MakeDb/"
    log:
        "logs/{imgt_release}/changeo_makedb/{sample}_{time}_{chain}.log"
    shell:
        '''
            set -e

            # log stdout and stderr
            mkdir -p logs/{rule}
            exec &> >(tee {log})
            
            mkdir -p {params.out_dir}

            MakeDb.py igblast \
                -i {input.blastn_fmt7} \
                -r {input.gapped_fasta} \
                -s {input.mod_clones} \
                --extended \
                --partial \
                --outdir {params.out_dir}
        '''





# rule complex_conversion:
#     input:
#         "{dataset}/inputfile"
#     output:
#         "{dataset}/file.{group}.txt"
#     shell:
#         "somecommand --group {wildcards.group} < {input} > {output}"


# rule processing_step:
#     input:
#         # [...]
#     output:
#         # [...]
#     run:
#         commands = [
#             "somecommand {input} > tempfile",
#             "othercommand tempfile {output}"
#         ]
#         for c in commands:
#             shell(c)

# rule aggregate:
#     input:
#         expand("{dataset}/a.{ext}", dataset=DATASETS, ext=FORMATS)
#     output:
#         "aggregated.txt"
#     shell:
#         ...
# If FORMATS=["txt", "csv"] 
