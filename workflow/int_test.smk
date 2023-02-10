####################################################################
###   Integration Testing for Extracellular RNA Bioinformatics   ###
####################################################################

import pandas as pd
import re
container: config["container"]

# Setup sample name index as a python dictionary

libraries = pd.read_table(config["data_dir"] + "/inputs/libraries.tsv")
library_indict = libraries["library"].tolist()
file_indict = libraries["file"].tolist()
lib_dict = dict(zip(library_indict, file_indict))

LIBRARIES = list(lib_dict.keys())
FASTQS = list(lib_dict.values())

rule all:
    input:
        #expand(config["data_dir"] + "/fastq/raw/{library}.fastq.gz", library = lib_dict.keys()),
        #expand(config["data_dir"] + "/fastq/trim/{library}.fastq.gz", library = LIBRARIES),
        #expand(config["data_dir"] + "/bam/{align_step}/{library}_{align_step}_ReadsPerGene.out.tab", library = LIBRARIES, align_step = ["mirna"]),
        #expand(config["data_dir"] + "/counts/{library}_{align_step}_counts.tsv", library = LIBRARIES, align_step = ["mirna"]),
        #config["data_dir"] + "/counts/counts.tsv",
        config["data_dir"] + "/de/de.Rdata",

rule symlink_inputs:
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        config["data_dir"] + "/fastq/raw/{library}.fastq.gz"
    shell:
        """
        ln -sf --relative {input} {output}
        """

include: "exrna.smk"
