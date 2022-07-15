#!/usr/bin/env python3
import os
import collections
from Bio import SeqIO
import csv
from itertools import groupby, count
# from cigar import Cigar

from piranha.utils.config import *

def id_reference_cns(aln):
    ref_seq = False
    for record in SeqIO.parse(aln, KEY_FASTA):
        if not ref_seq:
            ref_seq = str(record.seq).upper()
        else:
            cns_seq = str(record.seq).upper()
    
    return ref_seq,cns_seq

def merge_indels(indel_list,prefix):
    if indel_list:
        groups = groupby(indel_list, key=lambda item, c=count():item-next(c))
        tmp = [list(g) for k, g in groups]
        merged_indels = []
        for i in tmp:
            indel = f"{i[0]}:{prefix}{len(i)}"
            merged_indels.append(indel)
        return merged_indels
    
    return indel_list

def find_variants(reference_seq,query_seq):

    non_amb = ["A","T","G","C"]
    
    variants = []

    snps =[]
    insertions = []
    deletions = []

    for i in range(len(query_seq)):
        # bases[0] == query seq
        # bases[1] == ref seq
        bases = [query_seq[i],reference_seq[i]]
        if bases[0] != bases[1]:
            
            if bases[0] in non_amb and bases[1] in non_amb:
                #if neither of them are ambiguous
                snp = f"{i+1}:{bases[1]}{bases[0]}" # position-reference-query
                snps.append(snp)
            elif bases[0]=='-':
                #if there's a gap in the query, means a deletion
                deletions.append(i+1)
            elif bases[1]=='-':
                #if there's a gap in the ref, means an insertion
                insertions.append(i+1)

    insertions = merge_indels(insertions,"ins")
    deletions = merge_indels(deletions,"del")

    for var_list in [snps,insertions,deletions]:
        for var in var_list:
            variants.append(var)

    variants = sorted(variants, key = lambda x : int(x.split(":")[0]))

    return variants

def find_ambiguity_pcent(query_seq):
    ambiguity_count = 0
    for i in query_seq:
        if i not in ["A","T","G","C","-"]:
            ambiguity_count +=1
    ambiguity_pcent = round(100*(ambiguity_count/len(query_seq)),2)
    return ambiguity_pcent

def parse_variants(alignment,out_report,barcode,reference):

    ref_seq,cns_seq = id_reference_cns(alignment)
    variants = find_variants(ref_seq,cns_seq)
    var_string = ";".join(variants)
    with open(out_report,"w") as fw:
        fw.write(f"{barcode},{reference},{len(variants)},{var_string}\n")

def join_variant_files(header_fields,in_files,output):
    with open(output,"w") as fw:
        header = ",".join(header_fields) + "\n"
        fw.write(header)
        for in_file in in_files:
            with open(in_file, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    fw.write(f"{l}\n")

# def adjust_position(cs,primer_length):

#     # new cigar string
#     c = Cigar(cs)
#     items = list(c.items())
#     start_mask = items[0][0]
#     end_mask = items[-1][0]

#     amended_start_mask = start_mask + primer_length
#     amended_end_mask = end_mask + primer_length

#     new_cigar = c.mask_left(amended_start_mask).mask_right(amended_end_mask)

#     return new_cigar.cigar

# def soft_mask_primer_sites(input_sam, output_sam, primer_length):
#     with open(output_sam, "w") as fw:
#         with open(input_sam,"r") as f:
#             for l in f:
#                 l=l.rstrip("\n")

#                 if not l.startswith("@"):
#                     tokens = l.split("\t")
#                     new_tokens = tokens
#                     new_tokens[5] = adjust_position(tokens[5],primer_length)
#                     new_l = "\t".join(new_tokens)
#                     fw.write(f"{new_l}\n")
#                 else:
#                     fw.write(f"{l}\n")



def get_variation_pcent(ref,stats):
    ref_record = SeqIO.read(ref,KEY_FASTA)
    
    #set up site counter, 0-based
    variant_sites = collections.Counter()

    c = 0
    for record in SeqIO.parse(fasta,KEY_FASTA):
        c +=1
        index = 1
        for pos in zip(str(record.seq),str(ref_record.seq)):
            
            if pos[0] != pos[1]:
                variant_sites[index]+=1

            index +=1

    
    std_list = ["A","T","C","G","-"]
    variant_info = {}
    for site in variant_sites:
        variant_info[site] = {}
        for base in std_list:
            variant_info[site][base] = 0

    for record in SeqIO.parse(fasta,KEY_FASTA):
        for site in variant_sites:

            variant = str(record.seq)[site-1]
            if variant not in std_list:
                variant_info[site]["N"]+=1
            else:
                variant_info[site][variant]+=1

    x = []
    y = []
    info = []
    variation_info = []
    for site in variant_sites:
        
        pcent_variants = round(100*(variant_sites[site]/c), 1)
        x = site
        y = pcent_variants

        var_counts = variant_info[site]
        site_data = {"Position":x,"Percentage":y}
        for i in var_counts:
            site_data[f"{i} reads"] = var_counts[i]
        
        variation_info.append(site_data)

    return variation_info