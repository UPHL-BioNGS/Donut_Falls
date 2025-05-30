#!/usr/bin/env python3
import glob
import json
import csv
from os.path import exists
def busco_results():
    dict = {}
    files = glob.glob('short_summary*txt')
    for file in files:
      sample_analysis = file.split('.')[-2]
      with open(file, 'r') as f:
        for line in f:
          if 'C:' and 'S:' and 'D:' and 'F:' and 'M:' and 'n:' in line:
            dict[sample_analysis] = line.strip()
            break
    return dict

def circulocov_results():
    dict = {}
    files = glob.glob('*overall_summary.txt')
    for file in files:
      sample_analysis = file.replace('_overall_summary.txt', '').replace('_reoriented', '')
      dict[sample_analysis] = {}
      with open(file, 'r') as f:
        for line in f:
          parts = line.split()
          if parts[2] == 'all':
            dict[sample_analysis]['coverage'] = parts[7]
            dict[sample_analysis]['nanopore_numreads'] = parts[5]
            if len(parts) >= 9:
              dict[sample_analysis]['illumina_numreads'] = parts[9]
          if parts[2] == 'missing':
            if len(parts) > 8:
              unmapped_illumina = parts[8]
            else:
              unmapped_illumina = 0
            dict[sample_analysis]['unmapped_nanopore'] = parts[4]
            dict[sample_analysis]['unmapped_illumina'] = unmapped_illumina
    return dict

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
        key = row[header] + '_' + row[header2]
        dict[row[header]][key] = row
    return dict

def final_file(dict):
    with open('donut_falls_summary.json', 'w') as json_file:
      json.dump(dict, json_file, indent=4)

def mash_file(file):
    dict = {}
    with open(file, mode = 'r') as file:
      reader = csv.DictReader(file, delimiter='\\t', fieldnames=['sample', 'illumina','nanopore', 'dist', 'pvalue', 'hash'])
      for row in reader:
        key = row['sample']
        dict[key] = row
    return dict

def tsv_file(dict):
    final_dict = {}
    with open('donut_falls_summary.tsv', 'w') as tsv:
      i = 0
      sorted_keys = sorted(dict.keys())
      for key in sorted_keys:
        final_dict[key] = {}
        for result_key in ['name', 'number_of_reads', 'mean_read_length', 'mean_qual', 'total_illumina_reads', 'nanopore_illumina_mash_distance', 'assemblers']:
          final_dict[key][result_key] = dict[key][result_key]
        
        results = ['total_length', 'num_contigs', 'circ_contigs', 'coverage', 'unmapped_nanopore', 'unmapped_nanopore_pc', 'unmapped_illumina', 'unmapped_illumina_pc']
        for assembler in ['flye', 'raven']:
          if assembler in dict[key]['assemblers']:
            if assembler in dict[key].keys():
              for result in results + ['busco']:
                final_dict[key][assembler + '_' + result] = dict[key][assembler][result]
              final_dict[key][assembler + '_busco_polished'] = dict[key][assembler]['busco_pypolca']
              final_dict[key][assembler + '_quality_before_polishing'] = dict[key][assembler]['Consensus_Quality_Before_Polishing']
              final_dict[key][assembler + '_QV_before_polishing'] = dict[key][assembler]['Consensus_QV_Before_Polishing']
            else:
              for result in results + ['quality_before_polishing', 'QV_before_polishing' ]:
                final_dict[key][assembler + '_' + result] = 0
              for result in ['busco', 'busco_polished']:
                final_dict[key][assembler + '_' + result] = 'NF'
        if 'unicycler' in dict[key]['assemblers']:
          if 'unicycler' in dict[key].keys():
            for result in results + [ 'busco']:
              final_dict[key]['unicycler_' + result] = dict[key]['unicycler'][result]
          else:
            for result in results:
              final_dict[key]['unicycler_' + result] = 0
            final_dict[key]['unicycler_busco'] = 'NF'
        w = csv.DictWriter(tsv, final_dict[key].keys(), delimiter='\\t')
        if i < 1 :
          w.writeheader()
          i = i+1
        w.writerow(final_dict[key])

def main():

    seqkit_dict   = file_to_dict('file') if exists('file') else {}
    mash_dict     = mash_file('mash_summary.tsv') if exists('mash_summary.tsv') else {}
    pypolca_dict  = file_to_dict('pypolca_summary.tsv', 'sample', '\t') if exists('pypolca_summary.tsv') else {}
    gfastats_dict = file_to_dict_uniq('gfastats_summary.csv', 'sample', 'Header', ',') if exists('gfastats_summary.csv') else {}
    
    busco_dict = busco_results()
    circulocov_dict = circulocov_results()
    final_results = {}
    assemblers = ['flye', 'myloasm', 'raven', 'unicycler']

    for key in seqkit_dict.keys():
      if key.endswith('_illumina'):
        actual_key = key.replace('_illumina', '')
        if actual_key not in final_results.keys():
          final_results[actual_key] = {}
        final_results[actual_key]['total_illumina_reads'] = int(seqkit_dict[key]['number_of_reads'])
      else:
        if key not in final_results.keys():
          final_results[key] = {}
        if 'total_illumina_reads' not in final_results[key].keys():
          final_results[key]['total_illumina_reads'] = 0
        final_results[key]['name'] = key
        # from seqkit
        final_results[key]['number_of_reads']  = int(seqkit_dict[key]['number_of_reads'])
        final_results[key]['mean_read_length'] = float(seqkit_dict[key]['mean_read_length'])
        final_results[key]['mean_qual']        = float(seqkit_dict[key]['mean_qual'])
        # from mash
        if key in mash_dict.keys():
          final_results[key]['nanopore_illumina_mash_distance'] = float(mash_dict[key]['dist'])
        else:
          final_results[key]['nanopore_illumina_mash_distance'] = 'NF'
        final_results[key]['assemblers'] = '${params.assembler}'
        # for each assembler
        for assembler in assemblers:
          if key + '_' + assembler in gfastats_dict.keys():
            final_results[key][assembler] = {}
            final_results[key][assembler]['assembler'] = assembler
            # gfastats results
            total_length = 0
            num_circular = 0
            for contig in gfastats_dict[key + '_' + assembler].keys():
              total_length = total_length + int(gfastats_dict[key + '_' + assembler][contig]['Total segment length'])
              if gfastats_dict[key + '_' + assembler][contig]['Is circular'] == 'Y':
                num_circular = num_circular + 1
            final_results[key][assembler]['total_length'] = total_length
            final_results[key][assembler]['num_contigs']  = len(gfastats_dict[key + '_' + assembler].keys())
            final_results[key][assembler]['circ_contigs'] = num_circular
            # circulocov results
            if key + '_' + assembler in circulocov_dict.keys():
              if 'coverage' in circulocov_dict[key + '_' + assembler].keys():
                final_results[key][assembler]['coverage']          = float(circulocov_dict[key + '_' + assembler]['coverage'])
              else:
                final_results[key][assembler]['coverage']          = 'NF'
              if 'unmapped_nanopore' in circulocov_dict[key + '_' + assembler].keys():
                final_results[key][assembler]['unmapped_nanopore']    = int(circulocov_dict[key + '_' + assembler]['unmapped_nanopore'])
                final_results[key][assembler]['unmapped_nanopore_pc'] = float('{:.2f}'.format(int(final_results[key][assembler]['unmapped_nanopore']) / int(seqkit_dict[key]['number_of_reads']) * 100))
              else:
                final_results[key][assembler]['unmapped_nanopore']    = 'NF'
                final_results[key][assembler]['unmapped_nanopore_pc'] = 'NF'
              if 'unmapped_illumina' in circulocov_dict[key + '_' + assembler].keys():
                final_results[key][assembler]['unmapped_illumina'] = int(circulocov_dict[key + '_' + assembler]['unmapped_illumina'])
                if 'total_illumina_reads' in final_results[key].keys() and final_results[key]['total_illumina_reads'] > 0:
                  final_results[key][assembler]['unmapped_illumina_pc'] = float('{:.2f}'.format(int(final_results[key][assembler]['unmapped_illumina']) / int(final_results[key]['total_illumina_reads']) * 100 ))
                else:
                  final_results[key][assembler]['unmapped_illumina_pc'] = 0.0
              else:
                final_results[key][assembler]['unmapped_illumina'] = 'NF'
                final_results[key][assembler]['unmapped_illumina_pc'] = 'NF'
            # busco results
            if key + '_' + assembler in busco_dict.keys():
              final_results[key][assembler]['busco'] = busco_dict[key + '_' + assembler]
            elif key + '_' + assembler + '_reoriented' in busco_dict.keys():
              final_results[key][assembler]['busco'] = busco_dict[key + '_' + assembler + '_reoriented']
            else:
              final_results[key][assembler]['busco'] = 'NF'
            if assembler != 'unicycler':
              for step in ['polypolish', 'pypolca', 'medaka']:
                if key + '_' + assembler + '_' + step in busco_dict.keys():
                  final_results[key][assembler]['busco_' + step ] = busco_dict[key + '_' + assembler + '_' + step]
                else:
                  final_results[key][assembler]['busco_' + step ] = 'NF'
            # pypolca results
            if key + '_' + assembler in pypolca_dict.keys():
              if 'Consensus_Quality_Before_Polishing' in pypolca_dict[key + '_' + assembler].keys():
                final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = float(pypolca_dict[key + '_' + assembler]['Consensus_Quality_Before_Polishing'])
              else:
                final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = 'NF'
              if 'Consensus_QV_Before_Polishing' in pypolca_dict[key + '_' + assembler].keys():
                final_results[key][assembler]['Consensus_QV_Before_Polishing']      = float(pypolca_dict[key + '_' + assembler]['Consensus_QV_Before_Polishing'])
              else:
                final_results[key][assembler]['Consensus_QV_Before_Polishing']      = 'NF'
            elif assembler != 'unicycler':
              final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = 0.0
              final_results[key][assembler]['Consensus_QV_Before_Polishing']      = 0.0
    final_file(final_results)
    tsv_file(final_results)
if __name__ == '__main__':
  main()