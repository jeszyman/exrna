####################################################################
###   Integration Testing for Extracellular RNA Bioinformatics   ###
####################################################################

import pandas as pd
import re
container: config["container"]


# Setup sample name index as a python dictionary

libraries = pd.read_table(config["data_dir"] + "/inputs/libraries.tsv")
library_indict = libraries["library"].tolist()
file_indict = libraries["file_id"].tolist()
file_indict = [re.sub('$','.fastq.gz', file) for file in file_indict]
file_indict = [re.sub('^', config["data_dir"] + "/inputs/", file) for file in file_indict]
lib_dict = dict(zip(library_indict, file_indict))

LIBRARIES = list(lib_dict.keys())
FASTQS = list(lib_dict.values())

rule all:
    input:
        expand(config["data_dir"] + "/fastq/raw/{library}.fastq.gz", library = lib_dict.keys()),
        expand(config["data_dir"] + "/fastq/trim/{library}.fastq.gz", library = LIBRARIES),
        expand(config["data_dir"] + "/bam/{align_step}/{library}_{align_step}_ReadsPerGene.out.tab", library = LIBRARIES, align_step = ["mirna"]),
        expand(config["data_dir"] + "/counts/{library}_{align_step}_counts.tsv", library = LIBRARIES, align_step = ["mirna"]),
        config["data_dir"] + "/counts/counts.tsv",
        config["data_dir"] + "/de/de.Rdata",

include: "exrna.smk"
