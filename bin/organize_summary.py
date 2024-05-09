#!/usr/bin/env python3

import glob
import json
import csv
from os.path import exists

def busco_results():
    dict = {}
    files = glob.glob("short_summary*txt")
    for file in files:
      sample_analysis = file.split(".")[-2]
      with open(file, 'r') as f:
        for line in f:
          if "C:" and "S:" and "D:" and "F:" and "M:" and "n:" in line:
            dict[sample_analysis] = line.strip()
            break
    return dict

def circulocov_results():
    dict = {}
    files = glob.glob("*overall_summary.txt")
    for file in files:
      sample_analysis = file.replace("_overall_summary.txt", "").replace("_reoriented", "")
      dict[sample_analysis] = {}
      with open(file, 'r') as f:
        for line in f:
          parts = line.split()
          if parts[2] == "all":
            dict[sample_analysis]["coverage"] = parts[7]

          if parts[2] == "missing":
            if len(parts) > 8:
              unmapped_illumina = parts[8]
            else:
              unmapped_illumina = 0

            dict[sample_analysis]["unmapped_nanopore"] = parts[4]
            dict[sample_analysis]["unmapped_illumina"] = unmapped_illumina
    return dict

def fastp_results():
    dict = {}
    files = glob.glob("*_fastp_sr.json")
    for file in files:
      sample = file.replace('_fastp_sr.json', '')
      with open(file, 'r') as f:
        data = json.load(f)
        total_reads = data['summary']['before_filtering']['total_reads']
        dict[sample] = total_reads
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
        key = row[header] + "_" + row[header2]
        dict[row[header]][key] = row
    return dict

def final_file(dict):
    with open('donut_falls_summary.json', 'w') as json_file:
      json.dump(dict, json_file, indent=4)

def mash_file(file):
    dict = {}
    with open(file, mode = 'r') as file:
      reader = csv.DictReader(file, delimiter='\t', fieldnames=['illumina','nanopore', 'dist', 'pvalue', 'hash'])
      for row in reader:
        key = row['nanopore'].replace('.fastq.gz', '')
        dict[key] = row
    return dict

def tsv_file(dict):
    final_dict = {}
    with open('donut_falls_summary.tsv', 'w') as tsv:
      i = 0
      sorted_keys = sorted(dict.keys())
      for key in sorted_keys:
        final_dict[key] = {}
        final_dict[key]['name'] = dict[key]['name']
        final_dict[key]["number_of_reads"] = dict[key]['number_of_reads']
        final_dict[key]["mean_read_length"] = dict[key]['mean_read_length']
        final_dict[key]["mean_qual"] = dict[key]['mean_qual']
        final_dict[key]["total_illumina_reads"] = dict[key]['total_illumina_reads']
        final_dict[key]["nanopore_illumina_mash_distance"] = dict[key]['nanopore_illumina_mash_distance']
        final_dict[key]["assemblers"] = dict[key]['assemblers']
        
        if "flye" in dict[key]['assemblers'].replace("dragonflye","dragon"):
          if "flye" in dict[key].keys():
            final_dict[key]['flye_total_length'] = dict[key]['flye']['total_length']
            final_dict[key]['flye_num_contigs'] = dict[key]['flye']['num_contigs']
            final_dict[key]['flye_circ_contigs'] = dict[key]['flye']['circ_contigs']
            final_dict[key]['flye_coverage'] = dict[key]['flye']['coverage']
            final_dict[key]['flye_unmapped_nanopore'] = dict[key]['flye']['unmapped_nanopore']
            final_dict[key]['flye_unmapped_nanopore_pc'] = dict[key]['flye']['unmapped_nanopore_pc']
            final_dict[key]['flye_unmapped_illumina'] = dict[key]['flye']['unmapped_illumina']
            final_dict[key]['flye_unmapped_illumina_pc'] = dict[key]['flye']['unmapped_illumina_pc']
            final_dict[key]['flye_busco'] = dict[key]['flye']['busco']
            final_dict[key]['flye_busco_polished'] = dict[key]['flye']['busco_pypolca']
            final_dict[key]['flye_quality_before_polishing'] = dict[key]['flye']['Consensus_Quality_Before_Polishing']
            final_dict[key]['flye_QV_before_polishing'] = dict[key]['flye']['Consensus_QV_Before_Polishing']
          else:
            final_dict[key]['flye_total_length'] = 0
            final_dict[key]['flye_num_contigs'] = 0
            final_dict[key]['flye_circ_contigs'] = 0
            final_dict[key]['flye_coverage'] = 0
            final_dict[key]['flye_unmapped_nanopore'] = 0
            final_dict[key]['flye_unmapped_nanopore_pc'] = 0
            final_dict[key]['flye_unmapped_illumina'] = 0
            final_dict[key]['flye_unmapped_illumina_pc'] = 0
            final_dict[key]['flye_busco'] = "NF"
            final_dict[key]['flye_busco_polished'] = "NF"
            final_dict[key]['flye_quality_before_polishing'] = 0
            final_dict[key]['flye_QV_before_polishing'] = 0            



        if "raven" in dict[key]['assemblers']:
          if "raven" in dict[key].keys():
            final_dict[key]['raven_total_length'] = dict[key]['raven']['total_length']
            final_dict[key]['raven_num_contigs'] = dict[key]['raven']['num_contigs']
            final_dict[key]['raven_circ_contigs'] = dict[key]['raven']['circ_contigs']
            final_dict[key]['raven_coverage'] = dict[key]['raven']['coverage']
            final_dict[key]['raven_unmapped_nanopore'] = dict[key]['raven']['unmapped_nanopore']
            final_dict[key]['raven_unmapped_nanopore_pc'] = dict[key]['raven']['unmapped_nanopore_pc']
            final_dict[key]['raven_unmapped_illumina'] = dict[key]['raven']['unmapped_illumina']
            final_dict[key]['raven_unmapped_illumina_pc'] = dict[key]['raven']['unmapped_illumina_pc']
            final_dict[key]['raven_busco'] = dict[key]['raven']['busco']
            final_dict[key]['raven_busco_polished'] = dict[key]['raven']['busco_pypolca']
            final_dict[key]['raven_quality_before_polishing'] = dict[key]['raven']['Consensus_Quality_Before_Polishing']
            final_dict[key]['raven_QV_before_polishing'] = dict[key]['raven']['Consensus_QV_Before_Polishing']
          else:
            final_dict[key]['raven_total_length'] = 0
            final_dict[key]['raven_num_contigs'] = 0
            final_dict[key]['raven_circ_contigs'] = 0
            final_dict[key]['raven_coverage'] = 0
            final_dict[key]['raven_unmapped_nanopore'] = 0
            final_dict[key]['raven_unmapped_nanopore_pc'] = 0
            final_dict[key]['raven_unmapped_illumina'] = 0
            final_dict[key]['raven_unmapped_illumina_pc'] = 0
            final_dict[key]['raven_busco'] = "NF"
            final_dict[key]['raven_busco_polished'] = "NF"
            final_dict[key]['raven_quality_before_polishing'] = 0
            final_dict[key]['raven_QV_before_polishing'] = 0            

        if "unicycler" in dict[key]['assemblers']:
          if "unicycler" in dict[key].keys():
            final_dict[key]['unicycler_total_length'] = dict[key]['unicycler']['total_length']
            final_dict[key]['unicycler_num_contigs'] = dict[key]['unicycler']['num_contigs']
            final_dict[key]['unicycler_circ_contigs'] = dict[key]['unicycler']['circ_contigs']
            final_dict[key]['unicycler_coverage'] = dict[key]['unicycler']['coverage']
            final_dict[key]['unicycler_unmapped_nanopore'] = dict[key]['unicycler']['unmapped_nanopore']
            final_dict[key]['unicycler_unmapped_nanopore_pc'] = dict[key]['unicycler']['unmapped_nanopore_pc']
            final_dict[key]['unicycler_unmapped_illumina'] = dict[key]['unicycler']['unmapped_illumina']
            final_dict[key]['unicycler_unmapped_illumina_pc'] = dict[key]['unicycler']['unmapped_illumina_pc']
            final_dict[key]['unicycler_busco'] = dict[key]['unicycler']['busco']
          else:
            final_dict[key]['unicycler_total_length'] = 0
            final_dict[key]['unicycler_num_contigs'] = 0
            final_dict[key]['unicycler_circ_contigs'] = 0
            final_dict[key]['unicycler_coverage'] = 0
            final_dict[key]['unicycler_unmapped_nanopore'] = 0
            final_dict[key]['unicycler_unmapped_nanopore_pc'] = 0
            final_dict[key]['unicycler_unmapped_illumina'] = 0
            final_dict[key]['unicycler_unmapped_illumina_pc'] = 0
            final_dict[key]['unicycler_busco'] = "NF"

        w = csv.DictWriter(tsv, final_dict[key].keys(), delimiter='\t')
        if i < 1 :
          w.writeheader()
          i = i+1
        w.writerow(final_dict[key])

def main():
    if exists('nanoplot_summary.csv') :
      nanoplot_dict = file_to_dict('nanoplot_summary.csv', 'sample', ',')
    else:
      nanoplot_dict = {}

    if exists('mash_summary.tsv') :
      mash_dict = mash_file('mash_summary.tsv')
    else:
      mash_dict = {}

    if exists('pypolca_summary.tsv') :
      pypolca_dict  = file_to_dict('pypolca_summary.tsv', 'sample', '\t')
    else:
      pypolca_dict = {}

    if exists('gfastats_summary.csv') :
      gfastats_dict = file_to_dict_uniq('gfastats_summary.csv', 'sample', 'Header', ',')
    else:
      gfastats_dict = {}

    fastp_dict = fastp_results()

    busco_dict = busco_results()

    circulocov_dict = circulocov_results()

    final_results = {}
    assemblers = ['dragonflye', 'flye', 'hybracter', 'raven', 'unicycler']
    for key in nanoplot_dict.keys():
      final_results[key] = {}
      final_results[key]['name'] = key

      # from nanostas
      final_results[key]['number_of_reads']  = nanoplot_dict[key]['number_of_reads']
      final_results[key]['mean_read_length'] = nanoplot_dict[key]['mean_read_length']
      final_results[key]['mean_qual']        = nanoplot_dict[key]['mean_qual']

      # from fastp
      if key in fastp_dict.keys():
        final_results[key]['total_illumina_reads'] = fastp_dict[key]
      else:
        final_results[key]['total_illumina_reads'] = 0

      # from mash
      if key in mash_dict.keys():
        final_results[key]['nanopore_illumina_mash_distance'] = mash_dict[key]['dist']
      else:
        final_results[key]['nanopore_illumina_mash_distance'] = "NF"

      final_results[key]["assemblers"] = "flye,raven,unicycler"

      # for each assembler
      for assembler in assemblers:
        if key + "_" + assembler in gfastats_dict.keys():
          final_results[key][assembler] = {}
          final_results[key][assembler]['assembler'] = assembler

          # gfastats results
          total_length = 0
          num_circular = 0
          for contig in gfastats_dict[key + "_" + assembler].keys():
            total_length = total_length + int(gfastats_dict[key + "_" + assembler][contig]["Total segment length"])
            if gfastats_dict[key + "_" + assembler][contig]["circular"] == "Y":
              num_circular = num_circular + 1

          final_results[key][assembler]['total_length'] = total_length
          final_results[key][assembler]['num_contigs']  = len(gfastats_dict[key + "_" + assembler].keys())
          final_results[key][assembler]['circ_contigs'] = num_circular

          # circulocov results
          if key + "_" + assembler in circulocov_dict.keys():
            if 'coverage' in circulocov_dict[key + '_' + assembler].keys():
              final_results[key][assembler]['coverage']          = circulocov_dict[key + '_' + assembler]['coverage']
            else:
              final_results[key][assembler]['coverage']          = "NF"

            if 'unmapped_nanopore' in circulocov_dict[key + '_' + assembler].keys():
              final_results[key][assembler]['unmapped_nanopore']    = circulocov_dict[key + '_' + assembler]['unmapped_nanopore']
              final_results[key][assembler]['unmapped_nanopore_pc'] = "{:.2f}".format(int(final_results[key][assembler]['unmapped_nanopore']) / int(nanoplot_dict[key]['number_of_reads']) * 100)
            else:
              final_results[key][assembler]['unmapped_nanopore']    = "NF"
              final_results[key][assembler]['unmapped_nanopore_pc'] = "NF"

            if 'unmapped_illumina' in circulocov_dict[key + '_' + assembler].keys():
              final_results[key][assembler]['unmapped_illumina'] = circulocov_dict[key + '_' + assembler]['unmapped_illumina']
              if 'total_illumina_reads' in final_results[key].keys() and final_results[key]['total_illumina_reads'] > 0:
                final_results[key][assembler]['unmapped_illumina_pc'] = "{:.2f}".format(int(final_results[key][assembler]['unmapped_illumina']) / int(final_results[key]['total_illumina_reads']) * 100 )
              else:
                final_results[key][assembler]['unmapped_illumina_pc'] = 0.0
            else:
              final_results[key][assembler]['unmapped_illumina'] = "NF"
              final_results[key][assembler]['unmapped_illumina_pc'] = "NF"

          # busco results
          if key + "_" + assembler in busco_dict.keys():
            final_results[key][assembler]['busco'] = busco_dict[key + "_" + assembler]
          elif key + "_" + assembler + '_reoriented' in busco_dict.keys():
            final_results[key][assembler]['busco'] = busco_dict[key + "_" + assembler + '_reoriented']
          else:
            final_results[key][assembler]['busco'] = "NF"

          if assembler != 'unicycler':
            for step in ['polypolish', 'pypolca', 'medaka']:
              if key + "_" + assembler + '_' + step in busco_dict.keys():
                final_results[key][assembler]['busco_' + step ] = busco_dict[key + "_" + assembler + '_' + step]
              else:
                final_results[key][assembler]['busco_' + step ] = 'NF'

          # pypolca results
          if key + "_" + assembler in pypolca_dict.keys():
            if 'Consensus_Quality_Before_Polishing' in pypolca_dict[key + "_" + assembler].keys():
              final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = pypolca_dict[key + "_" + assembler]['Consensus_Quality_Before_Polishing']
            else:
              final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = "NF"
            if 'Consensus_QV_Before_Polishing' in pypolca_dict[key + "_" + assembler].keys():
              final_results[key][assembler]['Consensus_QV_Before_Polishing']      = pypolca_dict[key + "_" + assembler]['Consensus_QV_Before_Polishing']
            else:
              final_results[key][assembler]['Consensus_QV_Before_Polishing']      = "NF"

          elif assembler != 'unicycler':
            final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = 0
            final_results[key][assembler]['Consensus_QV_Before_Polishing']      = 0

    final_file(final_results)
    tsv_file(final_results)

if __name__ == "__main__":
  main()
