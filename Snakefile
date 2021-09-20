configfile: "config.yaml"

pipeline = "gbs-genetic-fitness" # replace with your pipeline's name


include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
READS = config["READS"]
PREFIX = config["PREFIX"]
GENOME_SIZE = config["GENOME_SIZE"]

idx_list = [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]

localrules: one_line_fasta, change_fq_extension
rule all:
    input:
        files_log,
        expand("DATA/{prefix}_oneline.fa", prefix=PREFIX),
        expand('DATA/{prefix}.fq.gz', prefix=PREFIX),
        expand('{prefix}.longstitch.done', prefix=PREFIX),
        expand("{prefix}.fa{idx}", prefix=PREFIX, idx = idx_list)


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
            shell("cp {input} {params.prefix}.fq.gz")

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
print(LONGSTITCH)
rule prep_ref:
    input:
        LONGSTITCH[0] + ".ntLink-arks.longstitch-scaffolds.fa"
    output:
        multiext("{prefix}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
    message:
        'Rule {rule} processing'
    shell:
        """
module load samtools
bwa index -a bwtsw {input}
samtools faidx {input}
        """