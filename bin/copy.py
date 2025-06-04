#!/usr/bin/env python3
import glob
import os
import shutil

def check_circ(file, name):
    for x in [".fasta", "_reoriented", "_clair3", "_polypolish", "_pypolca"]:
        name = name.replace(x, "")

    with open(file, "r") as f:
        for line in f:
            line = line.strip()
            searching_line = line.split(",")
            if searching_line[0] == name and searching_line[-1] == "N":
                return False
    return True


def sub_fasta(fasta):
    with open(fasta, 'r') as file:
        i = 0
        with open(f"consensus/sub_{fasta}", 'w') as outfile:
            for line in file:
                line = line.strip()
                if line.startswith('>') and i < 1:
                    outfile.write(f">{line.split()[0]} [location=chromosome][topology=circular][completeness=complete]\n")
                    i += 1
                elif line.startswith('>') and i >= 1:
                    outfile.write(f">{line.split()[0]} [plasmid-name=unnamed{i}][topology=circular][completeness=complete]\n")
                else:
                    outfile.write(f"{line}\n")


def main():
  
    os.mkdir('consensus')
  
    fasta = glob.glob('*.fasta')[0]
  
    shutil.copy(fasta, f"consensus/{fasta}")
  
    if check_circ('gfastats_summary.csv', fasta):
   
        sub_fasta(fasta)

if __name__ == '__main__':
    main()