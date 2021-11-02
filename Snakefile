configfile: "config_meagre.yaml"

pipeline = "gbs-genetic-fitness"

include: "rules/create_file_log.smk"

ASSEMBLY = config["ASSEMBLY"]
LONGREADS = config["LONGREADS"]
PREFIX = config["PREFIX"]
GENOME_SIZE = config["GENOME_SIZE"]
PARAMETERS_FILE = config["PARAMETERS_FILE"]
RUN_SCAFFOLDING = config["RUN_SCAFFOLDING"]

idx_list = [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]


# If COMPARISON_GENOME is set, compute whole genome alignment between polished assembly and comparison genome
if "COMPARISON_GENOME" in config:
    comp_genome_results = expand("results/genome_alignment/{prefix}_{species}.png", prefix=PREFIX, species = config["COMPARISON_GENOME"].keys())
    MIN_ALIGNMENT_LENGTH = config["MIN_ALIGNMENT_LENGTH"]
    MIN_QUERY_LENGTH = config["MIN_QUERY_LENGTH"]
else:
    comp_genome_results = []

if RUN_SCAFFOLDING == "N":
    extra_genome_stats = []
elif RUN_SCAFFOLDING == "Y":
    extra_genome_stats = f"results/assembly_stats_{PREFIX}_new.txt"


localrules: one_line_fasta, get_assembly_stats, plink_dist_matr, plink_IBD, bcftools_stats, plink_inbreeding, plot_pca, create_fastgbs_checkpoint_file, create_file_log, plot_aligned_genomes
rule all:
    input:
        files_log,
        expand("results/assembly_stats_{prefix}_{asm}.txt", asm = "ref", prefix=PREFIX),
        expand("results/{prefix}_FastGBS_platypus.recode.vcf.gz", prefix = PREFIX),
        expand("results/{prefix}_FastGBS_platypus.recode.vcf.gz.stats", prefix=PREFIX),
        expand("fitness/{prefix}.dist.gz", prefix=PREFIX),
        expand("fitness/{prefix}.genome.gz", prefix=PREFIX),
        expand("fitness/{prefix}.het.gz", prefix=PREFIX),
        expand('results/{prefix}_pca.html', prefix=PREFIX),
        expand("{prefix}_fast_gbs.done", prefix=PREFIX),
        comp_genome_results,
        extra_genome_stats
        

        # expand("{prefix}_fastgbs_check.done", prefix=PREFIX),
        # expand("results/{prefix}_FastGBS_platypus_log.txt", prefix = PREFIX),
        # expand("results/{prefix}.eigenvec", prefix=PREFIX),
        # expand("results/{prefix}.eigenval", prefix=PREFIX),
        # expand("{prefix}_oneline.fa", prefix=PREFIX),
        # expand("refgenome/{prefix}.fa{idx}", prefix=PREFIX, idx = idx_list),
        # expand("results/assembly_stats_{prefix}_{version}.txt", prefix = PREFIX, version = ["old", "new"]),


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

def run_scaffolding(RUN_SCAFFOLDING):
    if RUN_SCAFFOLDING == "Y":
        return("{prefix}_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa")
    elif RUN_SCAFFOLDING == "N":
        return(ASSEMBLY) 

rule prep_ref:
    input:
        run_scaffolding(RUN_SCAFFOLDING)
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

cp {input} {output.refgenome}
bwa index -a bwtsw {output.refgenome}
samtools faidx {output.refgenome}
        """

def assembly_stats_input(wildcards):
    if wildcards.asm == "ref":
        return(ASSEMBLY)
    if wildcards.asm == "new":
        return("refgenome/{prefix}.fa")


rule get_assembly_stats:
    input:
        assembly_stats_input
    output:
        "results/assembly_stats_{prefix}_{asm}.txt",
    message:
        'Rule {rule} processing'
    params:
        script = os.path.join(workflow.basedir, "scripts/get_assembly_stats.py"),
    log:
        err = "logs_slurm/get_assembly_stats_{prefix}_{asm}.err",
    shell:
        """
        python {params.script} {input} > {output} 2> {log.err}
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
printf 'IMPUTATION' > {output}
        """

rule run_fast_GBS2:
    input:
        parameters = PARAMETERS_FILE,
        refgenome = rules.prep_ref.output.refgenome,
        checkpoint = rules.create_fastgbs_checkpoint_file.output
    output:
        # vcf = "results/{prefix}_FastGBS_platypus.vcf",
        vcf_filtered ="results/{prefix}_FastGBS_platypus.recode.vcf",
        # log = "results/{prefix}_FastGBS_platypus_log.txt",
        check = touch("{prefix}_fast_gbs.done"),
        # GT = "results/{prefix}_FastGBS_platypus.recode.GT.FORMAT",
    message:
        'Rule {rule} processing'
    conda:
        "envs/fastgbs2_dependencies.yaml"
    shell:
        """
module load htslib/gcc/64/1.9 samtools bwa vcftools bedtools
export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/Platypus_0.8.1:$PATH
# export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/srg-extractor:$PATH

./fastgbs_V2.sh {input.parameters}
        """

rule zip_index_vcf:
    input:
        rules.run_fast_GBS2.output.check,
        vcf = "results/{prefix}_FastGBS_platypus.recode.vcf"
        # vcf = rules.run_fast_GBS2.output.vcf_filtered
    output:
        vcf = "results/{prefix}_FastGBS_platypus.recode.vcf.gz",
        idx = "results/{prefix}_FastGBS_platypus.recode.vcf.gz.tbi"
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
        # rules.run_fast_GBS2.output.check,
        vcf = rules.zip_index_vcf.output.vcf, 
        idx = rules.zip_index_vcf.output.idx
    output:
        "results/{prefix}_FastGBS_platypus.recode.vcf.gz.stats"
    message:
        'Rule {rule} processing'
    shell:
        """
        module load bcftools
        bcftools stats -s - {input.vcf} > {output}
        """


rule plink_dist_matr:
    input:
        vcf = rules.zip_index_vcf.output.vcf, 
        idx = rules.zip_index_vcf.output.idx
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
        vcf = rules.zip_index_vcf.output.vcf, 
        idx = rules.zip_index_vcf.output.idx
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
        vcf = rules.zip_index_vcf.output.vcf, 
        idx = rules.zip_index_vcf.output.idx
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
        vcf = rules.zip_index_vcf.output.vcf,
        idx = rules.zip_index_vcf.output.idx
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

        
def get_ref_path(wildcards):
    '''
    Get genome path for comparison species
    '''
    return(config["COMPARISON_GENOME"][wildcards.species])


rule align_genomes:
    input:
        assembly = rules.prep_ref.output.refgenome,
        comparison = get_ref_path
    output:
        "results/genome_alignment/{prefix}_vs_{species}.paf"
    message:
        'Rule {rule} processing'
    shell:
        """
minimap2 -t 12 -cx asm5 {input.comparison} {input.assembly} > {output}
        """

rule plot_aligned_genomes:
    input:
        rules.align_genomes.output
    output:
        "results/genome_alignment/{prefix}_{species}.png"
    message:
        'Rule {rule} processing'
    params:
        script = os.path.join(workflow.basedir, "scripts/pafCoordsDotPlotly.R"),
        min_alignment_length = MIN_ALIGNMENT_LENGTH,
        min_query_length = MIN_QUERY_LENGTH, 
        outdir = "results/genome_alignment/"
    shell:
        """
        module load R
        Rscript {params.script} -i {input} -o {wildcards.prefix}_{wildcards.species} -s -t -x -m {params.min_alignment_length} -q {params.min_query_length} -l
        mv {wildcards.prefix}_{wildcards.species}.png {params.outdir}
        """

onsuccess:
    print("Workflow done. Cleaning up....")
    shell("rm data/*.fastq")
