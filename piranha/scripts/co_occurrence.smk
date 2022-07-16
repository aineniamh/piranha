import os
import collections
from Bio import SeqIO
import yaml
import json

from piranha.analysis.get_co_occurrence import *

from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *



BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]


rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"co_occurrence_info.json")


rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")


rule get_co_occurrence:
    input:
        fasta = expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","pseudoaln.fasta"),reference=REFERENCES),
        csv = os.path.join(config[KEY_TEMPDIR],"variants.csv")
    output:
        json = os.path.join(config[KEY_TEMPDIR],"co_occurrence_info.json")
    run:
        counter = {}
        with open(input.csv,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                # if "Sabin" in row["reference"]:
                variants = row["variants"]
                ref_file = ""
                for ref_file in input.fasta:
                    if row["reference"] in ref_file:

                        read_fasta_file = ref_file
                if variants:
                    counter[row["reference"]] = get_combinations(variants,read_fasta_file,row["reference"],BARCODE,10)
        
        with open(output.json, "w") as fw:
            fw.write(json.dumps(counter))