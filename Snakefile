configfile: "config.yaml"

pipeline = "gbs-genetic-fitness" # replace with your pipeline's name


include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
READS = config["READS"]
PREFIX = config["PREFIX"]
GENOME_SIZE = config["GENOME_SIZE"]
PARAMETERS_FILE = config["PARAMETERS_FILE"]

idx_list = [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]


BARCODES, = glob_wildcards("fast-gbs_v2/barcodes/{barcodes}")

localrules: one_line_fasta, change_fq_extension
rule all:
    input:
        files_log,
        expand("DATA/{prefix}_oneline.fa", prefix=PREFIX),
        expand('DATA/{prefix}.fq.gz', prefix=PREFIX),
        expand('{prefix}.longstitch.done', prefix=PREFIX),
        expand("fast-gbs_v2/refgenome/{prefix}.fa{idx}", prefix=PREFIX, idx = idx_list),
        expand('fast-gbs_v2/barcodes/barcodes_{barcodes}', barcodes=BARCODES),
        expand("fast-gbs_v2/results/{prefix}_FastGBS_platypus.vcf", prefix=PREFIX)


rule one_line_fasta:
    input:
        ASSEMBLY
    output:
        "DATA/{prefix}_oneline.fa"
    message:
        'Rule {rule} processing'
    shell:
        'seqtk seq -l0 {input} > {output}'

rule change_fq_extension:
    input:
        READS
    output:
        'DATA/{prefix}.fq.gz'
    message:
        'Rule {rule} processing'
    params:
        prefix=PREFIX
    run:
        if not input[0].endswith("fq.gz"):
            shell("cp {input} DATA/{params.prefix}.fq.gz")

rule scaffolding_long_reads:
    input:
        draft = rules.one_line_fasta.output, 
        reads = rules.change_fq_extension.output
    output:
        done = '{prefix}.longstitch.done',
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
    message:
        'Rule {rule} processing'
    params:
        size=GENOME_SIZE,
        draft = "DATA/{prefix}_oneline",
        reads = "DATA/{prefix}"
    shell:
        """
        longstitch ntLink-arks draft={params.draft} reads={params.reads} G={params.size}
        touch {output.done}
        """

LONGSTITCH, = glob_wildcards("{longstitch}.ntLink-arks.longstitch-scaffolds.fa")

rule prep_ref:
    input:
        longstitch_done = rules.scaffolding_long_reads.output.done,
        ref = LONGSTITCH[0] + ".ntLink-arks.longstitch-scaffolds.fa"
    output:
        refgenome = "fast-gbs_v2/refgenome/{prefix}.fa",
        idx = multiext("fast-gbs_v2/refgenome/{prefix}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
    message:
        'Rule {rule} processing'
    shell:
        """
cp {input.ref} {output.refgenome}
module load samtools
bwa index -a bwtsw {output.refgenome}
samtools faidx {output.refgenome}
        """

rule convert_barcode_format:
    input:
        "fast-gbs_v2/barcodes/{barcodes}"
    output:
        'fast-gbs_v2/barcodes/barcodes_{barcodes}'
    message:
        'Rule {rule} processing'
    shell:
        """
sh fast-gbs_v2/txt2unix.sh {input} > {output}
rm {input}
        """

rule run_fast_GBS2:
    input:
        parameters = PARAMETERS_FILE,
        barcodes = expand('fast-gbs_v2/barcodes/barcodes_{barcodes}', barcodes=BARCODES),
        refgenome = rules.prep_ref.output
    output:
        # "fast-gbs_v2/results/List.bam",
        "fast-gbs_v2/results/{prefix}_FastGBS_platypus.vcf",
        "fast-gbs_v2/results/{prefix}_FastGBS_platypus.recode.vcf",
        "fast-gbs_v2/results/{prefix}_FastGBS_platypus.GT.FORMAT",
        "fast-gbs_v2/results/{prefix}_FastGBS_platypus_recode_imputed.vcf.gz"
    message:
        'Rule {rule} processing'
    shell:
        './fast-gbs_v2/fastgbs_V2.sh {input.parameters}'
