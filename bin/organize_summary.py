#!/usr/bin/env python3

import glob
import json
import csv
from os.path import exists

def file_to_dict(file, header, delim):
    dict = {}
    with open(file, mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter=delim)
        for row in reader:
            key = row[header]
            dict[key] = row
    return dict

def file_to_dict_uniq(file, header, header2, delim):
    dict = {}
    with open(file, mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter=delim)
        for row in reader:
            if row[header] not in dict.keys():
                dict[row[header]] = {}

            key = row[header] + "_" + row[header2]
            dict[row[header]][key] = row
    return dict


def final_file(dict, assemblers):
    with open('donut_falls_summary.json', 'w') as json_file:
        json.dump(dict, json_file, indent=4)

def main():
    
    if exists('nanoplot_summary.csv') :
        nanoplot_dict = file_to_dict('nanoplot_summary.csv', 'sample', ',')
    
    if exists('pypolca_summary.tsv') :
        pypolca_dict  = file_to_dict('pypolca_summary.tsv', 'sample', '\t')

    if exists('gfastats_summary.csv') :
        gfastats_dict = file_to_dict_uniq('gfastats_summary.csv', 'sample', 'Header', ',')

    busco_dict = {}
    busco_files = glob.glob("short_summary*txt")
    for file in busco_files:
        sample_analysis = file.split(".")[-2]
        with open(file, 'r') as f:
            for line in f:
                if "C:" and "S:" and "D:" and "F:" and "M:" and "n:" in line:
                    busco_dict[sample_analysis] = line.strip()
                    break

    circulocov_dict = {}
    circulocov_files = glob.glob("*overall_summary.txt")
    for file in circulocov_files:
        sample_analysis = file.replace("_overall_summary.txt", "").replace("_reoriented", "")
        circulocov_dict[sample_analysis] = {}
        with open(file, 'r') as f:
            for line in f:
                parts = line.split()
                if parts[2] == "all":
                    circulocov_dict[sample_analysis]["coverage"] = parts[7]

                if "missing" in line:
                    if len(parts) > 8:
                        unmapped_illumina = parts[8]
                    else:
                        unmapped_illumina = 0
        
                    circulocov_dict[sample_analysis]["unmapped_nanopore"] = parts[4]
                    circulocov_dict[sample_analysis]["unmapped_illumina"] = unmapped_illumina
    
    final_results = {}
    assemblers = ['dragonflye', 'flye', 'hybracter', 'raven', 'unicycler']
    for key in nanoplot_dict.keys():
        final_results[key] = {}
        final_results[key]['name'] = key

        # from nanostas
        final_results[key]['number_of_reads'] = nanoplot_dict[key]['number_of_reads']
        final_results[key]['mean_read_length'] = nanoplot_dict[key]['mean_read_length']
        final_results[key]['mean_qual'] = nanoplot_dict[key]['mean_qual']
        for assembler in assemblers:
            if key + "_" + assembler in gfastats_dict.keys():
                final_results[key][assembler] = {}

                # gfastats results
                total_length  = 0
                num_circular = 0
                for contig in gfastats_dict[key + "_" + assembler].keys():
                    total_length = total_length + int(gfastats_dict[key + "_" + assembler][contig]["Total segment length"])
                    if gfastats_dict[key + "_" + assembler][contig]["circular"] == "Y":
                        num_circular = num_circular + 1
                final_results[key][assembler]['total_length'] = total_length
                final_results[key][assembler]['num_contigs'] = len(gfastats_dict[key + "_" + assembler].keys())
                final_results[key][assembler]['circ_contigs'] = num_circular
                
                # circulocov results
                if key + "_" + assembler in circulocov_dict.keys():
                    final_results[key][assembler]['coverage'] = circulocov_dict[key + '_' + assembler]['coverage']
                    final_results[key][assembler]['unmapped_nanopore'] = circulocov_dict[key + '_' + assembler]['unmapped_nanopore']
                    final_results[key][assembler]['unmapped_illumina'] = circulocov_dict[key + '_' + assembler]['unmapped_illumina']

                # busco results
                if key + "_" + assembler in busco_dict.keys():
                    final_results[key][assembler]['busco'] = busco_dict[key + "_" + assembler]
                if key + "_" + assembler + '_reoriented' in busco_dict.keys():                
                    final_results[key][assembler]['busco'] = busco_dict[key + "_" + assembler + '_reoriented']
                for step in ['polypolish', 'pypolca', 'medaka']:
                    if key + "_" + assembler + '_' + step in busco_dict.keys():                
                        final_results[key][assembler]['busco_' + step ] = busco_dict[key + "_" + assembler + '_' + step]
                    else:
                        final_results[key][assembler]['busco_' + step ] = 'NF'

                # pypolca results
                if key + "_" + assembler in pypolca_dict.keys():
                    final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = pypolca_dict[key + "_" + assembler]['Consensus_Quality_Before_Polishing']
                    final_results[key][assembler]['Consensus_QV_Before_Polishing'] = pypolca_dict[key + "_" + assembler]['Consensus_QV_Before_Polishing']
                else:
                    final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = 0
                    final_results[key][assembler]['Consensus_QV_Before_Polishing'] = 0

    final_file(final_results, assemblers)

if __name__ == "__main__":
    main()