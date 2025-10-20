#!/usr/bin/env python3
import os
import shutil
import pandas as pd

def clean_name(sample_assembler):
    for x in [".fasta", "_reoriented", "_clair3", "_polypolish", "_pypolca"]:
        sample_assembler = sample_assembler.replace(x, '')

    sample, assembler = sample_assembler.rsplit('_', 1)
    return sample, assembler

def sub_fasta(fasta):
    with open(fasta, 'r') as file:
        i = 0
        with open(f"consensus/sub_{fasta}", 'w') as outfile:
            for line in file:
                line = line.strip()
                if line.startswith('>') and i < 1:
                    outfile.write(f"{line.split()[0]} [location=chromosome][topology=circular][completeness=complete]\n")
                    i += 1
                elif line.startswith('>') and i >= 1:
                    outfile.write(f"{line.split()[0]} [plasmid-name=unnamed{i}][topology=circular][completeness=complete]\n")
                else:
                    outfile.write(f"{line}\n")


#os.mkdir('consensus')

fasta = "test_unicycler.fasta"

shutil.copy(fasta, f"consensus/{fasta}")

sample, assembler = clean_name(fasta)
df = pd.read_table('assembly_info.csv')
df = df[(df['sample'] == sample) & (df['assembler'] == assembler)].copy()
if not (df['circ.'] == "N").any():
    sub_fasta(fasta)

