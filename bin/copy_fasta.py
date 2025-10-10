#!/usr/bin/env python3
import glob
import os
import pandas as pd

def clean_name(name):
    for x in [".fasta", "_reoriented", "_clair3", "_polypolish", "_pypolca"]:
        name = name.replace(x, '')
    return name

def copy_fasta(fasta, df):
    with open(fasta, 'r') as file:
        with open(f"consensus/{fasta}", 'w') as outfile:
            for line in file:
                line = line.strip().replace(".1,mult",",mult")
                if line.startswith('>'):
                    header = line.split()[0].replace(">", "")
                    length = df.loc[df['Header'] == header, 'Total segment length'].values[0]
                    circular = df.loc[df['Header'] == header, 'Is circular'].values[0]
                    gc_per = df.loc[df['Header'] == header, 'GC content %'].values[0]
                    outfile.write(f">{header} length={length} circular={circular} gc_per={gc_per}\n")
                else:
                    outfile.write(f"{line}\n")

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


def main():
    #os.mkdir('consensus')

    fasta = glob.glob('*.fasta')[0]
    name = clean_name(fasta)
    pd.set_option('display.max_colwidth', None)
    df = pd.read_csv('gfastats_summary.csv')
    df = df[df['sample'] == name].copy()
    df['Header'] = df['Header'].str.replace("duplicated-no.1","duplicated-no").str.replace("duplicated-yes.1","duplicated-yes").str.replace("duplicated-probably.1","duplicated-probably")
    df['Is circular'] = df['Is circular'].replace({"N": "false", "Y": "true"})

    copy_fasta(fasta, df)

    if not (df['Is circular'] == "N").any():
        sub_fasta(fasta)

if __name__ == '__main__':
    main()
