configfile: "config.yaml"

pipeline = "gbs-genetic-fitness" # replace with your pipeline's name


include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
LONGREADS = config["LONGREADS"]
PREFIX = config["PREFIX"]
GENOME_SIZE = config["GENOME_SIZE"]
PARAMETERS_FILE = config["PARAMETERS_FILE"]

idx_list = [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]


localrules: one_line_fasta, get_assembly_stats, plink_dist_matr, plink_IBD, bcftools_stats, plink_inbreeding, plot_pca
rule all:
    input:
        files_log,
        # expand("{prefix}_oneline.fa", prefix=PREFIX),
        # expand("refgenome/{prefix}.fa{idx}", prefix=PREFIX, idx = idx_list),
        # expand("results/assembly_stats_{prefix}_{version}.txt", prefix = PREFIX, version = ["old", "new"]),
        # expand("results/{prefix}_FastGBS_platypus.vcf", prefix = PREFIX),
        # expand("results/{prefix}_FastGBS_platypus_log.txt", prefix = PREFIX),
        expand("results/{prefix}_FastGBS_platypus.vcf.gz.stats", prefix=PREFIX),
        expand("fitness/{prefix}.dist.gz", prefix=PREFIX),
        expand("fitness/{prefix}.genome.gz", prefix=PREFIX),
        expand("{prefix}_fastgbs_check.done", prefix=PREFIX),
        expand("fitness/{prefix}.het.gz", prefix=PREFIX),
        expand("results/{prefix}.eigenvec", prefix=PREFIX),
        expand("results/{prefix}.eigenval", prefix=PREFIX),
        expand('results/{prefix}_pca.html', prefix=PREFIX)
        



rule one_line_fasta:
    input:
        ASSEMBLY
    output:
        "{prefix}_oneline.fa"
    message:
        'Rule {rule} processing'
    log:
        err = "logs_slurm/one_line_fasta_{prefix}.err",
    shell:
        'seqtk seq -l0 {input} > {output} 2> {log.err}'



rule scaffolding_long_reads:
    input:
        draft = rules.one_line_fasta.output,
        reads =  LONGREADS
    output:
        "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa"
    message:
        'Rule {rule} processing'
    params:
        size=GENOME_SIZE,
        draft = os.path.splitext(rules.one_line_fasta.output[0])[0],
        reads = os.path.splitext(os.path.splitext(os.path.basename(LONGREADS))[0])[0]
    shell:
        """
        longstitch ntLink-arks draft={params.draft} reads={params.reads} G={params.size}
        """


rule prep_ref:
    input:
        ref = "{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa"
    output:
        refgenome = "refgenome/{prefix}.fa",
        idx = multiext("refgenome/{prefix}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
    message:
        'Rule {rule} processing'
    shell:
        """
module load samtools
module load htslib
module load bwa

cp {input.ref} {output.refgenome}
bwa index -a bwtsw {output.refgenome}
samtools faidx {output.refgenome}
        """

rule get_assembly_stats:
    input:
        new = "refgenome/{prefix}.fa",
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
        err_new = "logs_slurm/get_assembly_stats_{prefix}_new.err",
    shell:
        """
        python {params.script} {input.new} > {output.new} 2> {log.err_new}
        python {params.script} {input.old} > {output.old} 2> {log.err_old}
        """


rule create_fastgbs_checkpoint_file:
    input:
        PARAMETERS_FILE
    output:
        "checkpoint_" + PARAMETERS_FILE
    message:
        'Rule {rule} processing'
    shell:
        """
printf 'IMPUTATION\nSUMMARY' > {output}
        """

rule run_fast_GBS2:
    input:
        parameters = PARAMETERS_FILE,
        refgenome = rules.prep_ref.output.refgenome,
        checkpoint = rules.create_fastgbs_checkpoint_file.output
    output:
        # vcf = "results/{prefix}_FastGBS_platypus.vcf",
        # "results/{prefix}_FastGBS_platypus.recode.vcf",
        # "results/{prefix}_FastGBS_platypus.GT.FORMAT",
        # log = "results/{prefix}_FastGBS_platypus_log.txt",
        check = touch("{prefix}_fast_gbs.done"),
    message:
        'Rule {rule} processing'
    conda:
        "envs/fastgbs2_dependencies.yaml"
    shell:
        """
        module load samtools bwa vcftools bedtools
        export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/Platypus_0.8.1:$PATH
        # export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/srg-extractor:$PATH

        ./fastgbs_V2.sh {input.parameters}
        """

rule index_vcf:
    input:
        rules.run_fast_GBS2.output.check,
        vcf = "results/{prefix}_FastGBS_platypus.vcf"
    output:
        vcf = "results/{prefix}_FastGBS_platypus.vcf.gz",
        idx = "results/{prefix}_FastGBS_platypus.vcf.gz.tbi"
    message:
        'Rule {rule} processing'
    shell:
        """
        module load samtools
        bgzip -c {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule bcftools_stats:
    input:
        rules.run_fast_GBS2.output.check,
        vcf = "results/{prefix}_FastGBS_platypus.vcf.gz",
        idx = "results/{prefix}_FastGBS_platypus.vcf.gz.tbi"
    output:
        "results/{prefix}_FastGBS_platypus.vcf.gz.stats"
    message:
        'Rule {rule} processing'
    shell:
        """
        module load bcftools
        bcftools stats -s - {input.vcf} > {output}
        """


rule plink_dist_matr:
    input:
        vcf = "results/{prefix}_FastGBS_platypus.vcf.gz",
        idx = "results/{prefix}_FastGBS_platypus.vcf.gz.tbi"
    output:
        "fitness/{prefix}.dist.gz"
    message:
        'Rule {rule} processing'
    params:
        prefix = os.path.join("fitness",PREFIX)
    shell:
        """
module load plink/1.9-180913        
plink --vcf {input.vcf} --distance gz --out {params.prefix} --double-id --allow-extra-chr
        """

rule plink_IBD:
    input:
        vcf = "results/{prefix}_FastGBS_platypus.vcf.gz",
        idx = "results/{prefix}_FastGBS_platypus.vcf.gz.tbi"
    output:
        "fitness/{prefix}.genome.gz"
    message:
        'Rule {rule} processing'
    params:
        prefix = os.path.join("fitness",PREFIX)
    shell:
        """
module load plink/1.9-180913 
plink --vcf {input.vcf} --genome "gz" --out {params.prefix} --double-id --allow-extra-chr
        """

rule plink_inbreeding:
    input:
        vcf = "results/{prefix}_FastGBS_platypus.vcf.gz",
        idx = "results/{prefix}_FastGBS_platypus.vcf.gz.tbi"
    output:
        "fitness/{prefix}.het.gz"
    message:
        'Rule {rule} processing'
    params:
        prefix = os.path.join("fitness",PREFIX)
    shell:
        """
module load plink/1.9-180913 
plink --vcf {input.vcf} --het "gz" --out {params.prefix} --double-id --allow-extra-chr
        """

rule PCA:
    input:
        vcf = "results/{prefix}_FastGBS_platypus.vcf.gz",
        idx = "results/{prefix}_FastGBS_platypus.vcf.gz.tbi"
    output:
        eigenvec = "results/{prefix}.eigenvec",
        eigenval = "results/{prefix}.eigenval",
        # pdf = "results/{prefix}.pdf"
    message:
        'Rule {rule} processing'
    params:
        prefix= os.path.join("results",PREFIX),
    shell:
        """
        module load plink/1.9-180913
        plink --vcf {input.vcf} --pca --double-id --out {params.prefix} --allow-extra-chr --threads 8
        """

rule plot_pca:
    input:
        eigenvec = "results/{prefix}.eigenvec",
        eigenval = "results/{prefix}.eigenval",
    output:
        'results/{prefix}_pca.html'
    message:
        'Rule {rule} processing'
    params:
        pyscript = os.path.join(workflow.basedir, "scripts/plot_interactive_pca.py"),
        output = os.path.join("results",PREFIX),
    shell:
        "python {params.pyscript} --eigenval {input.eigenval} --eigenvec {input.eigenvec} --out {params.output}"

        # Rscript {params.rscript} --eigenvec={output.eigenvec} --eigenval={output.eigenval} --output={output.pdf} --sample_list={input.sample_list}

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