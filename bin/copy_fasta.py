#!/usr/bin/env python3

import glob
import json
import csv
import os

def gfastats_to_dict(header_dict):
    dict = {}
    with open("gfastats_summary.csv", mode='r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row["sample"] == header_dict['name'] + "_" + header_dict['assembler']:
                key = row["Header"]

                dict[key] = row
    return dict

def circulocov_to_dict(header_dict):
    dict = {}
    with open("circulocov_summary.txt", mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:

            if row["sample"].replace("_reoriented","") == header_dict['name'] + "_" + header_dict['assembler'] :
                key = row["contigs"]

                dict[key] = row
    return dict

def copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict):
    with open(fasta, 'r') as file:
        with open(f"consensus/{header_dict['fasta']}", 'w') as outfile:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    contig = line.replace(">","").split()[0]
                    circular = gfa_dict[contig]['circular'].replace("Y","true").replace("N","false")
                    length = gfa_dict[contig]['Total segment length']
                    gc_per = gfa_dict[contig]['GC content %']
                    meandepth = circulocov_dict[contig]['nanopore_meandepth']
                    assembler = header_dict['assembler']
                    step = header_dict['step']
                    outfile.write(f">{contig} circ={circular} len={length} gc={gc_per} cov={meandepth} asmb={assembler} stp={step}\n")
                else:
                    outfile.write(f"{line}\n")

def main():

    os.mkdir("consensus")

    header_dict = {}
    fasta = glob.glob("*.fasta")[0]
    header_dict['fasta'] = fasta

    name = fasta.replace(".fasta", "")

    assemblers = ['dragonflye', 'flye', 'hybracter', 'raven', 'unicycler']
    steps = ["reoriented", 'polypolish', 'pypolca', 'medaka']
    for step in steps:
        if step in name:
            header_dict['step'] = step
            name = name.replace(f"_{step}","")
            break

    if 'step' not in header_dict.keys():
        header_dict['step'] = False

    for assembler in assemblers:
        if assembler in name:
            header_dict['assembler'] = assembler
            name = name.replace(f"_{assembler}","")
            break

    header_dict['name'] = name

    gfa_dict        = gfastats_to_dict(header_dict)
    circulocov_dict = circulocov_to_dict(header_dict)

    copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict)

if __name__ == "__main__":
    main()