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
        os.path.join(config[KEY_TEMPDIR],"consensus_sequences.fasta"),
        os.path.join(config[KEY_TEMPDIR],"variants.csv"),
        os.path.join(config[KEY_TEMPDIR],"masked_variants.csv"),
        expand(os.path.join(config[KEY_TEMPDIR],"snipit","{reference}.svg"), reference=REFERENCES),

rule files:
    params:
        ref=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.reference.fasta"),
        reads=os.path.join(config[KEY_TEMPDIR],"reference_groups","{reference}.fastq")

rule minimap2_racon:
    input:
        reads=rules.files.params.reads,
        ref=rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2_racon.log")
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.ref.sam")
    shell:
        """
        minimap2 -ax map-ont --score-N=0 --secondary=no {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """

rule racon:
    input:
        reads=rules.files.params.reads,
        fasta=rules.files.params.ref,
        sam= rules.minimap2_racon.output.sam
    output:
        os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","racon_cns.fasta")
    shell:
        "racon --no-trimming -t 1 {input.reads} {input.sam} {input.fasta} > {output}"

rule mafft_racon:
    input:
       fasta = rules.racon.output[0],
       ref = rules.files.params.ref
    output:
        temp_file = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","racon_ref.fasta"),
        aln = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","racon_ref.aln.fasta")
    shell:
        "cat {input.ref} {input.fasta} > {output.temp_file} && "
        "mafft {output.temp_file} > {output.aln} "

rule curate_indels_racon:
    input:
        aln = rules.mafft_racon.output.aln
    params:
        reference = "{reference}"
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","racon_cns.clean.fasta")
    run:
        clean_cns_gaps("racon", SAMPLE, input.aln, output.fasta)

rule minimap2_medaka:
    input:
        reads=rules.files.params.reads,
        ref=rules.curate_indels_racon.output.fasta
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2.log")
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.racon_cns.unmasked.sam")
    shell:
        """
        minimap2 -ax map-ont --score-N=0 --secondary=no {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """

# rule soft_mask_primers:
#     input:
#         sam = rules.minimap2_medaka.output.sam
#     output:
#         sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.racon_cns.masked.sam")
#     run:
#         # if config[KEY_ANALYSIS_MODE] == VALUE_ANALYSIS_MODE_WG_2TILE:
#         #     soft_mask_primer_sites(input.sam, output.sam,30)
#         # else:
#             # soft_mask_primer_sites(input.sam, output.sam, 30)
#         shell("cp {input.sam:q} {output.sam:q}")

rule sort_index:
    input:
        sam = rules.minimap2_medaka.output.sam
    output:
        bam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.sorted.bam"),
        index = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.sorted.bam.bai")
    shell:
        """
        samtools view -bS -F 4 {input.sam:q} | samtools sort -o {output[0]:q} &&
        samtools index {output.bam:q} {output.index:q}
        """

rule medaka_consensus:
    input:
        basecalls=rules.files.params.reads,
        draft=rules.curate_indels_racon.output.fasta,
        sam= rules.minimap2_medaka.output.sam,
        bam=rules.sort_index.output.bam
    params:
        outdir=os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka"),
        model = config[KEY_MEDAKA_MODEL]
    output:
        cns_mod = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","cns.mod.fasta"),
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","consensus_probs.hdf"),
        consensus= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","consensus.fasta")
    threads:
        workflow.cores
    shell:
        """
        [ ! -d {params.outdir:q} ] && mkdir {params.outdir:q}

        if [ -s {input.sam:q} ]

        then
            sed "s/[:,-]/_/g" {input.draft:q} > {output.cns_mod:q}
            medaka consensus --model "{params.model}" {input.bam:q} {output.probs:q} 
            medaka stitch {output.probs:q} {output.cns_mod:q} {output.consensus:q} 
        else
            touch {output.consensus:q}
            touch {output.probs:q}
            touch {output.cns_mod:q}
        fi
        """

"""
            medaka variant {input.draft:q} {output.probs:q} {output.vcf:q} 
            bgzip -f {output.vcf:q}
	        tabix -p vcf {output.vcf_gz:q}
            bcftools consensus {output.vcf_gz:q} {output.consensus:q}
"""


rule minimap2_cns:
    input:
        reads=rules.files.params.reads,
        ref=rules.medaka_consensus.output.consensus
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.minimap2_cns.log")
    output:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.consensus.unmasked.sam")
    shell:
        """
        minimap2 -ax map-ont --score-N=0 --secondary=no {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """

rule sort_index_cns:
    input:
        sam = rules.minimap2_cns.output.sam
    output:
        bam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped_cns.sorted.bam"),
        index = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped_cns.sorted.bam.bai")
    shell:
        """
        samtools view -bS -F 4 {input.sam:q} | samtools sort -o {output[0]:q} &&
        samtools index {output.bam:q} {output.index:q}
        """

rule medaka_cns_consensus:
    input:
        basecalls=rules.files.params.reads,
        draft=rules.medaka_consensus.output.consensus,
        sam= rules.minimap2_cns.output.sam,
        bam=rules.sort_index_cns.output.bam
    params:
        outdir=os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_cns"),
        model = config[KEY_MEDAKA_MODEL],
        reference = "{reference}"
    output:
        cns_mod = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_cns","cns.mod.fasta"),
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_cns","consensus_probs.hdf"),
        consensus= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_cns","consensus.fasta")
    threads:
        workflow.cores
    run:
        if "Sabin" not in params.reference:
            shell("""
                [ ! -d {params.outdir:q} ] && mkdir {params.outdir:q}

                if [ -s {input.sam:q} ]

                then
                    sed "s/[:,-]/_/g" {input.draft:q} > {output.cns_mod:q}
                    medaka consensus --model "{params.model}" {input.bam:q} {output.probs:q} 
                    medaka stitch {output.probs:q} {output.cns_mod:q} {output.consensus:q} 
                else
                    touch {output.consensus:q}
                    touch {output.probs:q}
                    touch {output.cns_mod:q}
                fi
                """)
        else:
            shell("""
                    touch {output.consensus:q}
                    touch {output.probs:q}
                    touch {output.cns_mod:q}
                """)


rule join_cns_ref:
    input:
        ref=rules.files.params.ref,
        medaka_cns=rules.medaka_consensus.output.consensus,
        cns_cns=rules.medaka_cns_consensus.output.consensus
    params:
        reference = "{reference}"
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_cns.fasta")
    run:
        with open(output[0],"w") as fw:
            display_name = ""
            ref_record = ""
            for record in SeqIO.parse(input.ref,KEY_FASTA):
                
                for field in record.description.split(" "):
                    if field.startswith("display_name"):
                        display_name = field.split("=")[1]
                ref_record = record

            if "Sabin" in params.reference:
                fw.write(f">{display_name} {ref_record.description}\n{ref_record.seq}\n")

                for record in SeqIO.parse(input.medaka_cns,KEY_FASTA):
                    record_name = str(SAMPLE).replace(" ","_")
                    fw.write(f">{record_name}\n{record.seq}\n")
            else:
                for record in SeqIO.parse(input.medaka_cns,KEY_FASTA):
                    fw.write(f">Consensus_{display_name}\n{record.seq}\n")

                for record in SeqIO.parse(input.cns_cns,KEY_FASTA):
                    record_name = str(SAMPLE).replace(" ","_")
                    fw.write(f">{record_name}\n{record.seq}\n")


rule align_cns_ref:
    input:
        rules.join_cns_ref.output.fasta
    output:
        aln = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_cns.aln.fasta")
    shell:
        """
        mafft {input:q} > {output:q}
        """

rule curate_variants:
    input:
        aln = rules.align_cns_ref.output.aln
    params:
        reference = "{reference}"
    output:
        masked = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","masked.csv"),
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka_cns_clean.fasta")
    run:
        masked = clean_medaka_cns(SAMPLE, input.aln, output.fasta)
        with open(output.masked,"w") as fw:
            for var in masked:
                site = int(var) + 1
                fw.write(f"{params.reference},{site},{masked[var]}\n")

rule gather_masked_variants:
    input:
        expand(rules.curate_variants.output.masked, reference=REFERENCES)
    output:
        masked = os.path.join(config[KEY_TEMPDIR],"masked_variants.csv")
    run:
        with open(output.masked,"w") as fw:
            fw.write("reference,site,variant\n")
            for i in input:
                with open(i, "r") as f:
                    for l in f:
                        l=l.rstrip("\n")
                        fw.write(l + '\n')


rule join_clean_cns_ref:
    input:
        ref=rules.files.params.ref,
        cns=rules.curate_variants.output.fasta,
        medaka_cns=rules.medaka_consensus.output.consensus
    params:
        reference = "{reference}"
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_medaka_cns_clean.fasta")
    run:
        with open(output[0],"w") as fw:
            display_name = ""
            ref_record = ""
            for record in SeqIO.parse(input.ref,KEY_FASTA):
                ref_record = record
                for field in record.description.split(" "):
                    if field.startswith("display_name"):
                        display_name = field.split("=")[1]

            if "Sabin" in params.reference:
                fw.write(f">{display_name} {ref_record.description}\n{ref_record.seq}\n")
            else:
                for record in SeqIO.parse(input.medaka_cns,KEY_FASTA):
                    fw.write(f">Consensus_{display_name} {record.description}\n{record.seq}\n")
            for record in SeqIO.parse(input.cns,KEY_FASTA):
                record_name = SAMPLE.replace(" ","_")
                fw.write(f">{record_name}\n{record.seq}\n")

rule align_clean_cns_ref:
    input:
        rules.join_clean_cns_ref.output.fasta
    output:
        aln = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_medaka_cns_clean.aln.fasta")
    shell:
        """
        mafft {input:q} > {output:q}
        """

rule make_snipit_graph:
    input:
        aln = rules.align_clean_cns_ref.output.aln
    params:
        out_stem = os.path.join(config[KEY_TEMPDIR],"snipit","{reference}")
    output:
        os.path.join(config[KEY_TEMPDIR],"snipit","{reference}.svg")
    run:
        try:
            shell("""snipit {input.aln:q} -o {params.out_stem} -f svg -c wes""")
        except:
            shell("touch {output[0]:q}")

rule assess_variants:
    input:
        rules.align_clean_cns_ref.output.aln
    params:
        reference = "{reference}"
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","variants.csv")
    run:
        parse_variants(input[0],output.csv,BARCODE,params.reference)

rule gather_variants:
    input:
        expand(rules.assess_variants.output.csv, reference=REFERENCES)
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"variants.csv")
    run:
        join_variant_files(VARIANT_CALLS_HEADER_FIELDS,input,output.csv)

rule gather_cns:
    input:
        variants = rules.gather_variants.output.csv,
        seqs = expand(rules.curate_variants.output.fasta, reference=REFERENCES)
    output:
        os.path.join(config[KEY_TEMPDIR],"consensus_sequences.fasta")
    run:
        var_dict = {}
        with open(input.variants, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                var_dict[row[KEY_REFERENCE]] = row

        with open(output[0],"w") as fw:
            for in_file in input.seqs:
                reference = os.path.basename(os.path.dirname(in_file))
                for record in SeqIO.parse(in_file,KEY_FASTA):
                    variant_count = var_dict[reference]["variant_count"]
                    variant_string = var_dict[reference]["variants"]
                    fw.write(f">{reference}|{BARCODE}|{variant_count}|{variant_string}\n{record.seq}\n")
