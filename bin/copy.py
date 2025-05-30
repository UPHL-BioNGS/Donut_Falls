#!/usr/bin/env python3
import glob
import csv
import os

def gfastats_to_dict(header_dict):
    result = {}
    with open('gfastats_summary.csv', mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row['sample'] == f"{header_dict['name']}_{header_dict['assembler']}":
                key = row['Header']
                result[key] = row
    return result

def circulocov_to_dict(header_dict):
    result = {}
    with open('circulocov_summary.txt', mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            if row['sample'].replace('_reoriented','') == f"{header_dict['name']}_{header_dict['assembler']}" :
                key = row['contigs']
                result[key] = row
    return result

def copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict):
    with open(fasta, 'r') as file:
        fasta_dict = {}

        circ_check = any(gfa_dict[contig].get('Is circular') == 'N' for contig in gfa_dict)

        for line in file:
            line = line.strip()
            if line.startswith('>'):
                contig    = str(line.replace('>','').split()[0])
                circular  = gfa_dict[contig]['Is circular'].replace('Y','true').replace('N','false')
                length    = gfa_dict[contig]['Total segment length']
                gc_per    = gfa_dict[contig]['GC content %']
                meandepth = circulocov_dict[contig]['nanopore_meandepth']
                assembler = header_dict['assembler']
                step      = header_dict['step']
                
                # creating the header
                header = f">{contig} circ={circular} len={length} gc={gc_per} cov={meandepth} asmb={assembler}"
                if assembler != 'unicycler':
                    header += f" stp={step}"
                header += "\n"
                # creating the dict
                fasta_dict[contig] = {
                    'seq' : '',
                    'header' : header,
                    'length' : int(length)
                }
            else:
                fasta_dict[contig]['seq'] += line
    
    sorted_dict = dict(sorted(fasta_dict.items(), key=lambda x: x[1]['length'], reverse = True))
    with open(f"consensus/{header_dict['fasta']}", 'w') as outfile:    
        for contig in sorted_dict:
            seq = '\n'.join([fasta_dict[contig]['seq'][i:i+70] for i in range(0, len(fasta_dict[contig]['seq']), 70)])
            outfile.write(fasta_dict[contig]['header'])
            outfile.write(f"seq\n")

    if not circ_check:
        with open(f"consensus/sub_{header_dict['fasta']}", 'w') as outfile:
            i = 0
            for contig in sorted_dict:
                if i < 1:
                    sub_header = f">{fasta_dict[contig]['header'].split()[0]} [location=chromosome][topology=circular][completeness=complete]\n"
                else:
                    sub_header = f">{fasta_dict[contig]['header'].split()[0]} [plasmid-name=unnamed{i}][topology=circular][completeness=complete]\n"

                i += 1

                outfile.write(sub_header)
                seq = '\n'.join([fasta_dict[contig]['seq'][i:i+70] for i in range(0, len(fasta_dict[contig]['seq']), 70)])
                outfile.write(seq + '\n')

def main():
  
  os.mkdir('consensus')
  header_dict = {}
  fasta = glob.glob('*.fasta')[0]
  header_dict['fasta'] = fasta
  name = fasta.replace('.fasta', '')
  assemblers = ['flye', 'myloasm', 'raven', 'unicycler']
  steps = ['reoriented', 'clair3', 'polypolish', 'pypolca']
  
  for step in steps:
    if step in name:
      header_dict['step'] = step
      name = name.replace(f"_{step}",'')
      break
  
    if 'step' not in header_dict.keys():
        header_dict['step'] = False

  for assembler in assemblers:
    if assembler in name:
      header_dict['assembler'] = assembler
      name = name.replace(f"_{assembler}",'')
      break

  header_dict['name'] = name
  gfa_dict        = gfastats_to_dict(header_dict)
  circulocov_dict = circulocov_to_dict(header_dict)
  copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict)

if __name__ == '__main__':
  main()