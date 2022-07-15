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
        cns=os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_medaka_cns_clean.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")

rule map_cns:
    input:
        cns = rules.files.params.cns,
        ref = rules.files.params.ref,
        reads = rules.files.params.reads
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2_cns.log")
    params:
        ref = "{reference}"
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"variation_analysis","{reference}","cns_mapped.sam")
    run:
        if params.ref.startswith("Sabin"):
            reference = input.ref
        else:
            reference = input.cns

        shell(f"minimap2 -ax map-ont \
                --score-N=0 \
                --secondary=no \
                -o '{output.sam}' \
                '{reference}' \
                '{input.reads}' \
                &> {log:q}")

rule sort_index:
    input:
        sam = rules.map_cns.output.sam
    output:
        bam = os.path.join(config[KEY_TEMPDIR],"variation_analysis","{reference}","cns_mapped.sorted.bam"),
        index = os.path.join(config[KEY_TEMPDIR],"variation_analysis","{reference}","cns_mapped.sorted.bam.bai")
    shell:
        """
        samtools view -bS -F 4 {input.sam:q} | samtools sort -o {output[0]:q} &&
        samtools index {output.bam:q} {output.index:q}
        """


rule pysamstats:
    input:
        cns = rules.files.params.cns,
        ref = rules.files.params.ref,
        bam = rules.sort_index.output.bam
    output:
        stats = os.path.join(config[KEY_TEMPDIR],"variation_analysis","{reference}","pysamstats.variation.tsv")
    run:
        if params.ref.startswith("Sabin"):
            reference = input.ref
        else:
            reference = input.cns
        shell(f"pysamstats -f '{reference}' -t '{input.bam}' > '{output.stats}'")

# rule sam_to_seq:
#     input:
#         sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.ref.sam"),
#         ref = rules.files.params.ref
#     log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
#     output:
#         fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","pseudoaln.fasta")
#     shell:
#         """
#         gofasta sam toMultiAlign -r {input.ref:q} -s {input.sam:q} -o {output[0]:q} &> {log}
#         """

rule sam_to_indels:
    input:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.ref.sam"),
        ref = rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
    output:
        ins = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","insertions.tsv"),
        dels = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","deletions.tsv")
    shell:
        """
        gofasta sam indels --threshold {config[min_read_depth]} -s {input.sam:q} --insertions-out {output.ins:q} --deletions-out {output.dels:q} 
        """


rule get_variation_info:
    input:
        expand(rules.files.params.reads, reference=REFERENCES),
        expand(rules.pysamstats.output.stats,reference=REFERENCES)
    output:
        json = os.path.join(config[KEY_TEMPDIR],"variation_info.json")
    run:
        # this is for making a figure
        variation_dict = {}
        for reference in REFERENCES:
            # ref = os.path.join(config[KEY_TEMPDIR],"reference_groups",f"{reference}.reference.fasta")
            stats = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{reference}","pysamstats.variation.tsv")
            
            var_dict = get_variation_pcent(ref,stats)
            variation_dict[reference] = var_dict

        with open(output.json, "w") as fw:
            fw.write(json.dumps(variation_dict))