"""
A snakemake pipeline to align and filter nucmer alignments of a fasta file query (1 or multiple sequences) on a fasta file subject (1 or multiple sequences).
For instance, a set of scaffolds (query) aligned to a set of chromosomes (subject).
"""


from snakemake.utils import min_version
from glob import glob

############################
## Minimal Snakemake version
############################
min_version("5.2.0")

###############
# Configuration
###############

configfile: "config.yaml"

TEMP_DIR = config["tempdir"]
RESULT_DIR  = config["resultdir"]

QUERIES_DIR = config["queries_dir"]
REFS_DIR = config["refs_dir"]

LENGTH_THRESHOLD = config["length_threshold"]
IDENTITY_THRESHOLD = config["identity_threshold"]
#N_THRESHOLD_ABS = config["n_threshold_abs"]
N_THRESHOLD_PERC = config["n_threshold_perc"]

###########
# Wildcards
###########
def get_queries(wildcards):
    "get wildcard values and returns a list of query fasta files"
    queries = glob.glob(QUERIES_DIR + "{wildcards.query}.fasta")
    return queries

def get_refs(wildcards):
    "get wildcard values and returns a list of reference fasta files"
    refs = glob.glob(REFS_DIR + "{wildcards.ref}.fasta")
    return refs

def get_list_of_fasta_file_names(directory,extension=".fasta"):
    "extracts the fasta file name based on a directory and a file extension"
    names = [f.split(extension)[0] for f in os.listdir(directory) if f.endswith(extension)]
    return names

QUERIES = get_list_of_fasta_file_names(QUERIES_DIR)
REFS = get_list_of_fasta_file_names(REFS_DIR)

###################################
# Docker container for the pipeline
###################################
# singularity: "docker://continuumio/miniconda3:4.4.10"
# needs to be fixed

##################
# Pipeline outputs
##################
rule all:
    input:
#        expand(RESULT_DIR + "{query}_vs_{ref}.sorted.coords",query=QUERIES,ref=REFS),
        RESULT_DIR + "results.tsv"
    message:"all done"


############
# Rules
#############

rule nucmer:
    input:
        query = QUERIES_DIR + "{query}.fasta",
        ref = REFS_DIR + "{ref}.fasta"
    output:
        RESULT_DIR + "delta/{query}_vs_{ref}.delta"
    message:
        "Starting nucmer alignment with {wildcards.query} on {wildcards.ref}"
    conda:
        "envs/mummer.yaml"
#    threads: 20
    params:
        query = "{query}",
        ref = "{ref}",
        prefix = RESULT_DIR + "delta/{query}_vs_{ref}"
    shell:
        "nucmer --threads={threads} --prefix={params.prefix} {input.ref} {input.query}"

rule get_N_locations:
    input:
        QUERIES_DIR + "{query}.fasta"
    output:
        temp(TEMP_DIR + "query_N/{query}.txt")
#    params:
#        n_threshold = N_THRESHOLD
    message:
        "Starting N location acquirement for {wildcards.query}"
    conda:
        "envs/percentage.yaml"
    shell:
        "Rscript "
        "scripts/get_N_locations.r "
        "--fasta {input} "
        "--out {output} "
#        "--n_threshold {params.n_threshold} "

rule delta_to_coords:
    input:
        RESULT_DIR + "delta/{query}_vs_{ref}.delta"
    output:
        temp(TEMP_DIR + "coords/{query}_vs_{ref}.coords")
    message:
        "converting {input} format to the coords format"
    conda:
        "envs/mummer.yaml"
#    threads: 20
    shell:
        "show-coords -r {input} > {output}"

rule sort_coords:
    input:
        TEMP_DIR + "coords/{query}_vs_{ref}.coords"
    output:
        temp(TEMP_DIR + "coords/{query}_vs_{ref}.sorted.coords")
    message:
        "sorting {input} file by coordinates"
#    threads: 20
    shell:
        "sort  -n -k4 {input} > {output}"

rule calculate_alignment_percentage:
    input:
#        lambda wildcards: [RESULT_DIR + "coords/{wildcards.query}_vs_{wildcards.ref}.sorted.coords"]
        coords = TEMP_DIR + "coords/{query}_vs_{ref}.sorted.coords",
        fasta = QUERIES_DIR + "{query}.fasta",
        nfile = TEMP_DIR + "query_N/{query}.txt"
    output:
        temp(TEMP_DIR + "{query}_vs_{ref}.txt")
    params:
        length_threshold = LENGTH_THRESHOLD,
        identity_threshold = IDENTITY_THRESHOLD,
#        n_threshold_abs = N_THRESHOLD_ABS,
        n_threshold_perc = N_THRESHOLD_PERC
    message:
        "calculating the percentage of aligned bases for {input}"
    conda:
        "envs/percentage.yaml"
    shell:
        "Rscript scripts/aligned_perc_calc.r "
        "--filename {input.coords} "
        "--fasta {input.fasta} "
        "--nfile {input.nfile} "
        "--out {output} "
        "--length_threshold {params.length_threshold} "
        "--identity_threshold {params.identity_threshold} "
#        "--n_threshold_abs {params.n_threshold_abs} "
        "--n_threshold_perc {params.n_threshold_perc} "

rule create_results_matrix:
    input:
        expand(TEMP_DIR + "{query}_vs_{ref}.txt", query=QUERIES, ref=REFS)
    output:
        RESULT_DIR + "results.tsv"
    message:
        "creating final results.tsv file"
    conda:
        "envs/merge.yaml"
    shell:
        "Rscript "
        "scripts/merge2matrix.r "
        "--filename {input} "
        "--out {output} "

# rule filter_alignments:
