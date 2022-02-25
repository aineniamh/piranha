import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *
from piranha.analysis.consensus_functions import *
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *


BARCODE = config[KEY_BARCODE]
SAMPLE = config[KEY_SAMPLE]
REFERENCES = config[BARCODE]

rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"consensus_sequences.fasta"),
        os.path.join(config[KEY_TEMPDIR],"variants.csv"),
        os.path.join(config[KEY_TEMPDIR],"variation_info.json"),
         os.path.join(config[KEY_TEMPDIR],"masked_variants.csv"),
        expand(os.path.join(config[KEY_TEMPDIR],"snipit","{reference}.svg"), reference=REFERENCES),
        expand(os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","deletions.tsv"), reference=REFERENCES)

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
        minimap2 -ax map-ont {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
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
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.racon_cns.sam")
    shell:
        """
        minimap2 -ax map-ont {input.ref:q} {input.reads:q} -o {output.sam:q} &> {log:q}
        """

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
        sam= rules.minimap2_medaka.output.sam
    params:
        outdir=os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka")
    output:
        cns_mod = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","cns.mod.fasta"),
        probs = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","consensus_probs.hdf"),
        consensus= os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","medaka","consensus.fasta")
    threads:
        2
    shell:
        """
        if [ -s {input.sam:q} ]
        
        then
            sed "s/[:,-]/_/g" {input.draft:q} > {output.cns_mod:q}
            medaka_consensus -i {input.basecalls:q} -d {output.cns_mod:q} -o {params.outdir:q} -t 2
        else
            touch {output.consensus:q}
            touch {output.probs:q}
            touch {output.cns_mod:q}
        fi
        """


rule sam_to_seq:
    input:
        sam = rules.minimap2_racon.output.sam,
        ref = rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","pseudoaln.fasta")
    shell:
        """
        gofasta sam toMultiAlign -r {input.ref:q} -s {input.sam:q} -o {output[0]:q} &> {log}
        """

rule sam_to_indels:
    input:
        sam = rules.minimap2_racon.output.sam,
        ref = rules.files.params.ref
    log: os.path.join(config[KEY_TEMPDIR],"logs","{reference}.gofasta.log")
    output:
        ins = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","insertions.tsv"),
        dels = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","deletions.tsv")
    shell:
        """
        gofasta sam indels --threshold {config[min_read_depth]} -s {input.sam:q} --insertions-out {output.ins:q} --deletions-out {output.dels:q} 
        """



rule join_cns_ref:
    input:
        ref=rules.files.params.ref,
        cns=rules.medaka_consensus.output.consensus
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_cns.fasta")
    run:
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input.ref,"fasta"):
                display_name = ""
                for field in record.description.split(" "):
                    if field.startswith("display_name"):
                        display_name = field.split("=")[1]

                fw.write(f">{display_name} {record.description}\n{record.seq}\n")
            for record in SeqIO.parse(input.cns,"fasta"):
                record_name = SAMPLE.replace(" ","_")
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

rule align_clean_cns_to_reference:
    input:
        ref = rules.files.params.ref,
        cns = rules.curate_variants.output.fasta
    params:
        sam = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","mapped.final.sam")
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","aligned_clean_cns.fasta")
    shell:
        """
        minimap2 -a -x asm20 --sam-hit-only --secondary=no -t  {workflow.cores} {input.ref:q} {input.cns} -o {params.sam:q} 
        gofasta sam toMultiAlign \
            -s {params.sam:q} \
            -t {workflow.cores} \
            --reference {input.ref:q} \
            > '{output.fasta}'
        """

rule get_variation_info:
    input:
        expand(rules.align_clean_cns_to_reference.output.fasta, reference=REFERENCES),
        expand(rules.files.params.reads, reference=REFERENCES),
        expand(rules.sam_to_seq.output.fasta,reference=REFERENCES)
    output:
        json = os.path.join(config[KEY_TEMPDIR],"variation_info.json")
    run:
        # this is for making a figure
        variation_dict = {}
        for reference in REFERENCES:
            ref = os.path.join(config[KEY_TEMPDIR],"reference_groups",f"{reference}.reference.fasta")
            fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{reference}","pseudoaln.fasta")
            cns = os.path.join(config[KEY_TEMPDIR],"reference_analysis",f"{reference}","aligned_clean_cns.fasta")
            
            var_dict = get_variation_pcent(ref,cns,fasta)
            variation_dict[reference] = var_dict

        with open(output.json, "w") as fw:
            fw.write(json.dumps(variation_dict))

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
        cns=rules.curate_variants.output.fasta
    output:
        fasta = os.path.join(config[KEY_TEMPDIR],"reference_analysis","{reference}","ref_medaka_cns_clean.fasta")
    run:
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input.ref,"fasta"):
                display_name = ""
                for field in record.description.split(" "):
                    if field.startswith("display_name"):
                        display_name = field.split("=")[1]

                fw.write(f">{display_name} {record.description}\n{record.seq}\n")
            for record in SeqIO.parse(input.cns,"fasta"):
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
                for record in SeqIO.parse(in_file,"fasta"):
                    variant_count = var_dict[reference]["variant_count"]
                    variant_string = var_dict[reference]["variants"]
                    fw.write(f">{reference}|{BARCODE}|{variant_count}|{variant_string}\n{record.seq}\n")