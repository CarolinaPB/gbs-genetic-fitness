configfile: "config.yaml"

pipeline = "gbs-genetic-fitness" # replace with your pipeline's name


include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
READS = config["READS"]
PREFIX = config["PREFIX"]
GENOME_SIZE = config["GENOME_SIZE"]
PARAMETERS_FILE = config["PARAMETERS_FILE"]

idx_list = [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]


localrules: one_line_fasta, get_assembly_stats
rule all:
    input:
        files_log,
        expand("DATA/{prefix}_oneline.fa", prefix=PREFIX),
        expand('{prefix}.longstitch.done', prefix=PREFIX),
        expand("refgenome/{prefix}.fa{idx}", prefix=PREFIX, idx = idx_list),
        expand("results/assembly_stats_{prefix}_{version}.txt", prefix = PREFIX, version = ["old", "new"]),
        expand("fastgbs_{prefix}.done", prefix=PREFIX),


rule one_line_fasta:
    input:
        ASSEMBLY
    output:
        "DATA/{prefix}_oneline.fa"
    message:
        'Rule {rule} processing'
    log:
        err = "logs_slurm/one_line_fasta_{prefix}.err",
        out = "logs_slurm/one_line_fasta_{prefix}.out"
    shell:
        'seqtk seq -l0 {input} > {output} 2> {log.err} 1> {log.out}'



rule scaffolding_long_reads:
    input:
        draft = rules.one_line_fasta.output,
        reads =  READS
    output:
        done = '{prefix}.longstitch.done',
    message:
        'Rule {rule} processing'
    params:
        size=GENOME_SIZE,
        draft = os.path.splitext(rules.one_line_fasta.output[0])[0],
        reads = os.path.splitext(os.path.splitext(READS)[0])[0]
    shell:
        """
        longstitch ntLink-arks draft={params.draft} reads={params.reads} G={params.size}
        touch {output.done}
        """

LONGSTITCH, = glob_wildcards("DATA/{longstitch}.ntLink-arks.longstitch-scaffolds.fa")

rule prep_ref:
    input:
        longstitch_done = rules.scaffolding_long_reads.output.done,
        ref = expand("DATA/{longstitch}.ntLink-arks.longstitch-scaffolds.fa", longstitch = LONGSTITCH)
    output:
        refgenome = "refgenome/{prefix}.fa",
        idx = multiext("refgenome/{prefix}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
    message:
        'Rule {rule} processing'
    shell:
        """
module load samtools
module load htslib

cp {input.ref} {output.refgenome}
bwa index -a bwtsw {output.refgenome}
samtools faidx {output.refgenome}
        """

rule get_assembly_stats:
    input:
        longstitch_done = rules.scaffolding_long_reads.output.done,
        new = expand("DATA/{longstitch}.ntLink-arks.longstitch-scaffolds.fa", longstitch = LONGSTITCH),
        old = ASSEMBLY
    output:
        old = "results/assembly_stats_{prefix}_old.txt",
        new = "results/assembly_stats_{prefix}_new.txt"
    message:
        'Rule {rule} processing'
    params:
        script = os.path.join(workflow.basedir, "scripts/get_assembly_stats.py"),
    log:
        err_old = "logs_slurm/get_assembly_stats_{prefix}_old.err",
        out_old = "logs_slurm/get_assembly_stats{prefix}_old.out",
        err_new = "logs_slurm/get_assembly_stats_{prefix}_new.err",
        out_new = "logs_slurm/get_assembly_stats{prefix}_new.out",
    shell:
        """
        python {params.script} {input.new} > {output.new} 2> {log.err_new} 1> {log.out_new}
        python {params.script} {input.old} > {output.old} 2> {log.err_old} 1> {log.out_old}
        """

rule run_fast_GBS2:
    input:
        parameters = PARAMETERS_FILE,
        # barcodes = expand('barcodes/barcodes_{barcodes}', barcodes=BARCODES),
        refgenome = rules.prep_ref.output.refgenome
    output:
        # "results/List.bam",
        # "results/{prefix}_FastGBS_platypus.vcf",
        # "results/{prefix}_FastGBS_platypus.recode.vcf",
        # "results/{prefix}_FastGBS_platypus.GT.FORMAT",
        # "results/{prefix}_FastGBS_platypus_recode_imputed.vcf.gz"
        "fastgbs_{prefix}.done"
    message:
        'Rule {rule} processing'
    conda:
        "envs/fastgbs_dependencies.yaml"
    shell:
        """
        module load samtools bwa vcftools htslib
        export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/Platypus_0.8.1:$PATH

        ./fastgbs_V2.sh {input.parameters} && touch {output}
        """



# rule cleanup_dirs:
#     input:
#         rules.run_fast_GBS2.output
#     output:
#         'cleanup_{prefix}.done'
#     message:
#         'Rule {rule} processing'
#     shell:
#         """
#         rm data/*.fastq
#         """

#to install platypus
# download from here
# module load python  python/2.7.15
# module load htslib

#https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data
# Place where you want to install
# tar -xvzf Platypus_x.x.x.tgz
# cd Platypus_x.x.x
# ./buildPlatypus.sh

# export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/Platypus_0.8.1:$PATH


# move fq files to data/ dir and rename to be 
# Paired-end reads:

#     FLOWCELL_LANE#_R1.fq.gz
#     FLOWCELL_LANE#_R2.fq.gz
# Single-end reads:

#     FLOWCELL_LANE#.fq.gz


# rule convert_barcode_format:
#     input:
#         "DATA/barcodes/{barcodes}"
#     output:
#         'barcodes/barcodes_{barcodes}'
#     message:
#         'Rule {rule} processing'
#     shell:
#         """
# sh txt2unix.sh {input} > {output}
        # """



        # "{prefix}_oneline.fa.{ls0}.tsv",
        # multiext("{prefix}_oneline.fa.{ls1}", ".ntLink.scaffolds.fa", ".stitch.abyss-scaffold.fa", ".stitch.abyss-scaffold.renamed.fa"),
        # multiext("{prefix}_oneline.fa.{longstitch}", 
        # ".dist.gv", 
        # "_main.tsv", 
        # "_original.gv", 
        # ".tigpair_checkpoint.tsv"),
        # multiext("{prefix}_oneline.fa.{longstitch}_{longstitch_2}", 
        # ".assembly_correspondence.tsv"
        # ".gv", 
        # ".log", 
        # ".scaffolds", 
        # ".scaffolds.fa"),
        # ls ="{ls2}.ntLink-arks.longstitch-scaffolds.fa"

        # rule change_fq_extension:
#     input:
#         READS
#     output:
#         'DATA/{prefix}.fq.gz'
#     message:
#         'Rule {rule} processing'
#     params:
#         prefix=PREFIX
#     run:
#         if not input[0].endswith("fq.gz"):
#             shell("cp {input} data/{params.prefix}.fq.gz")
#         else:
#             shell("cp {input} data/")