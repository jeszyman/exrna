############################
###    exRNA Pipeline    ###
############################

rule symlink_inputs:
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        config["data_dir"] + "/fastq/raw/{library}.fastq.gz"
    shell:
        """
        ln -sf --relative {input} {output}
        """

rule read_preprocessing:
    input:
        config["data_dir"] + "/fastq/raw/{library}.fastq.gz",
    params:
        script = config["exrna_script_dir"] + "/read_preprocessing.sh",
        threads = config["threads"]
    output:
        config["data_dir"] + "/fastq/trim/{library}.fastq.gz",
    resources:
        mem_mb=5000
    log:
        config["data_dir"] + "/logs/{library}_read_preprocessing.log",
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} \
        {log}
        """

rule make_star_mirna_index:
    input:
        fasta = config["data_dir"] + "/inputs/chr19.fa",
        gtf = config["data_dir"] + "/inputs/ENCFF470CZH.gtf",
    params:
        outdir = config["data_dir"] + "/ref/mirna_star",
        script = config["exrna_script_dir"] + "/make_star_mirna_index.sh",
        threads = config["threads"],
    output:
        done = touch(directory(config["data_dir"] + "/ref/mirna_star")),
    log:
        config["data_dir"] + "/logs/make_star_mirna_index.log",
    shell:
        """
        {params.script} \
        {input.fasta} \
        {input.gtf} \
        {params.outdir} \
        {params.threads} > {log}
        """

rule align_mirna:
    input:
        index = config["data_dir"] + "/ref/mirna_star",
        fq = config["data_dir"] + "/fastq/trim/{library}.fastq.gz",
    params:
        script = config["exrna_script_dir"] + "/align_mirna.sh",
        threads = config["threads"],
    output:
        config["data_dir"] + "/bam/mirna/{library}_mirna_ReadsPerGene.out.tab",
    log:
        config["data_dir"] + "/logs/{library}_align_mirna.log",
    shell:
        """
        {params.script} \
        {input.index} \
        {input.fq} \
        {params.threads} \
        {output} &> {log}
        """

rule transform_star_counts:
    input:
        config["data_dir"] + "/bam/{align_step}/{library}_{align_step}_ReadsPerGene.out.tab",
    params:
        script = config["exrna_script_dir"] + "/transform_star_counts.Rscript",
    output:
        config["data_dir"] + "/counts/{library}_{align_step}_counts.tsv",
    log:
        config["data_dir"] + "/logs/{library}_{align_step}_transform_star_counts.log",
    shell:
        """
        Rscript {params.script} \
	{input} \
	{output} \
        >& {log}
        """

rule merge_star_counts:
    input:
        expand(config["data_dir"] + "/counts/{library}_{align_step}_counts.tsv", library = LIBRARIES, align_step = ["mirna"]),
    params:
        script = config["exrna_script_dir"] + "/merge_star_counts.R",
    output:
        config["data_dir"] + "/counts/counts.tsv",
    log:
        config["data_dir"] + "/logs/merge_star_counts.log",
    shell:
        """
        Rscript {params.script} \
        "{input}" \
        {output} \
        > {log} 2>&1
        """

rule diff_express:
    input:
        counts = config["data_dir"] + "/counts/counts.tsv",
        coldata = config["data_dir"] + "/inputs/libraries.tsv",
    params:
        script = config["exrna_script_dir"] + "/diff_express.R",
        design = config["deseq_design"],
    output:
        config["data_dir"] + "/de/de.Rdata"
    log:
        config["data_dir"] + "/logs/diff_express.log",
    shell:
        """
        Rscript {params.script} \
        {input.counts} \
        {input.coldata} \
        "{params.design}" \
        {output} \
        > {log} 2>&1
        """
