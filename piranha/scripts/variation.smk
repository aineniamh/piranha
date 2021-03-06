import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *

from piranha.analysis.consensus_functions import *
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *


BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]


rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"variation_info.json"),
        expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","deletions.tsv"), reference=REFERENCES)


rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")


rule sam_to_seq:
    input:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.ref.sam"),
        sam_cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.consensus.unmasked.sam"),
        ref = rules.files.params.ref,
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","consensus.fasta")
    params:
        reference = "{reference}"
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","pseudoaln.fasta")
    run:
        if "Sabin" in params.reference:
            shell("""
                gofasta sam toMultiAlign -r {input.ref:q} -s {input.sam:q} -o {output[0]:q} &> {log}
                """)
        else:
            shell("""
                gofasta sam toMultiAlign -r {input.cns:q} -s {input.sam_cns:q} -o {output[0]:q} &> {log}
                """)


rule sam_to_indels:
    input:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.ref.sam"),
        sam_cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.consensus.unmasked.sam"),
        ref = rules.files.params.ref,
        cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","consensus.fasta")
    params:
        reference = "{reference}"
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
    output:
        ins = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","insertions.tsv"),
        dels = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","deletions.tsv")
    run:
        if "Sabin" in params.reference:
            shell("""
                gofasta sam indels --threshold {config[min_read_depth]} -s {input.sam:q} --insertions-out {output.ins:q} --deletions-out {output.dels:q} 
                """)
        else:
            shell("""
                gofasta sam indels --threshold {config[min_read_depth]} -s {input.sam_cns:q} --insertions-out {output.ins:q} --deletions-out {output.dels:q} 
                """)


rule get_variation_info:
    input:
        expand(rules.files.params.reads, reference=REFERENCES),
        expand(rules.sam_to_seq.output.fasta,reference=REFERENCES)
    output:
        json = os.path.join(config[KEY_TEMPDIR],"variation_info.json")
    run:
        # this is for making a figure
        variation_dict = {}
        for reference in REFERENCES:
            if "Sabin" in reference:
                ref = os.path.join(config[KEY_TEMPDIR],"reference_groups",f"{reference}.reference.fasta")
            else:
                ref = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{reference}","medaka","consensus.fasta")
            
            fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{reference}","pseudoaln.fasta")
            
            var_dict = get_variation_pcent(ref,fasta)
            variation_dict[reference] = var_dict

        with open(output.json, "w") as fw:
            fw.write(json.dumps(variation_dict))