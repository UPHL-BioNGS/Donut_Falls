#!/usr/bin/env python3
import glob
import json
import csv
from os.path import exists


def busco_results():
    results = {}
    files = glob.glob('short_summary*txt')
    for file in files:
      sample_assembler_step = file.split('.')[-2].split("_")
      if sample_assembler_step[-1] == "unicycler":
        assembler = "unicycler"
        step = "unicycler"
        sample = "_".join(sample_assembler_step[:-1])
      else:
        assembler = sample_assembler_step[-2]
        step = sample_assembler_step[-1]
        sample = "_".join(sample_assembler_step[:-2])
      
      if sample not in results.keys():
        results[sample] = {}
      if assembler not in results[sample].keys():
        results[sample][assembler] = { }      
      with open(file, 'r') as f:
        for line in f:
          if 'C:' and 'S:' and 'D:' and 'F:' and 'M:' and 'n:' in line:
            results[sample][assembler][step] = line.strip()
            break
    return results

def circulocov_results():
    results = {}
    files = glob.glob('*overall_summary.txt')
    for file in files:
      sample, assembler = file.replace('_reoriented', '').replace('_overall_summary.txt','').rsplit("_", 1)
      if sample not in results.keys():
        results[sample] = {}
      if assembler not in results[sample].keys():
        results[sample][assembler] = {}

      with open(file, 'r') as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
          results[sample][assembler][row['contigs']] = row

      results[sample][assembler]['all']['warnings'] = ""

      results[sample][assembler]['all']['unmapped_nanopore'] = results[sample][assembler]['missing']['nanopore_numreads'] if 'nanopore_numreads' in results[sample][assembler]['missing'].keys() else 0
      results[sample][assembler]['all']['unmapped_nanopore_pc'] = round(float(results[sample][assembler]['all']['unmapped_nanopore']) / float(results[sample][assembler]['all']['nanopore_numreads']), 2)
      if results[sample][assembler]['all']['unmapped_nanopore_pc'] > 0.1:
        results[sample][assembler]['all']['warnings'] += "High proportion of unmapped Nanopore reads,"

      if 'illumina_numreads' in results[sample][assembler]['missing'].keys():
        results[sample][assembler]['all']['unmapped_illumina'] = results[sample][assembler]['missing']['illumina_numreads']
        results[sample][assembler]['all']['unmapped_illumina_pc'] = round(float(results[sample][assembler]['all']['unmapped_illumina']) / float(results[sample][assembler]['all']['illumina_numreads']), 2)
        if results[sample][assembler]['all']['unmapped_illumina_pc'] > 0.1:
          results[sample][assembler]['all']['warnings'] += "High proportion of unmapped Illumina reads,"
    return results

def combine_results(seqkit_dict, mash_dict, pypolca_dict, assembly_info_dict, busco_dict, circulocov_dict):

    final_results = seqkit_dict.copy()

    tool_results = [
        ('mash', mash_dict),
        ('pypolca', pypolca_dict),
        ('assembly_info', assembly_info_dict),
        ('busco', busco_dict),
        ('circulocov', circulocov_dict),
    ]

    #final_results[key]['assemblers'] = '${params.assembler}'
    assemblers = 'flye,unicycler,raven,myloasm'

    for key in final_results:
        final_results[key]['assemblers'] = assemblers
        
        for tool_name, tool_dict in tool_results:
            if key in tool_dict:
                final_results[key][tool_name] = tool_dict[key]
    return final_results
        
def final_file(results):
    with open('donut_falls_summary.json', 'w') as json_file:
      json.dump(results, json_file, indent=4)

def assembly_info_file(file):
    results = {}
    with open(file, mode='r', newline='') as file:
      reader = csv.DictReader(file, delimiter="\t")
      for row in reader:
        sample = row['sample']
        assembler = row['assembler']
        if sample not in results.keys():
          results[sample] = { }
        if assembler not in results[sample].keys():
          results[sample][assembler] = {}
        results[sample][assembler][row['seq_name']] = row

    for sample in results.keys():
      for assembler in results[sample].keys():
        total_length = 0
        num_circular = 0
        contigs = list(results[sample][assembler].keys())
        results[sample][assembler]['num_contigs'] = len(contigs)
        for contig in contigs:
          total_length = total_length + int(results[sample][assembler][contig]['length'])
          if results[sample][assembler][contig]['circ.'] == 'Y':
            num_circular = num_circular + 1
        results[sample][assembler]['total_length'] = total_length
        results[sample][assembler]['circ_contigs'] = num_circular
    return results

def mash_file(file):
    results = {}
    with open(file, mode = 'r') as file:
      reader = csv.DictReader(file, delimiter='\t', fieldnames=['sample', 'illumina','nanopore', 'dist', 'pvalue', 'hash'])
      for row in reader:
        key = row['sample']
        results[key] = row
    return results

def pypolca_file(file):
    results = {}
    with open(file, mode='r', newline='') as file:
      reader = csv.DictReader(file, delimiter='\t')
      for row in reader:
        sample, assembler = row['sample'].rsplit("_", 1)
        results[sample] = { assembler : row }
    return results

def seqkit_file(file):
    results = {}
    with open(file, mode='r', newline='') as file:
      reader = csv.DictReader(file, delimiter='\t')
      for row in reader:
        key = row['sample']
        fastq = row['file']
        if key not in results.keys():
          results[key] = {}
        results[key][fastq] = row

    for key in results.keys():
      sorted_results = dict(sorted(results[key].items(), key=lambda fastq: float(fastq[1]['avg_len']), reverse=True))
      results[key] = sorted_results
      nanopore_fastq = list(results[key].keys())[0]
      results[key][nanopore_fastq]['platform'] = 'nanopore'
      results[key]['nanopore'] = results[key][nanopore_fastq]
      del results[key][nanopore_fastq]

      ill_fastq = []
      results_keys = list(results[key].keys())      
      for fastq in results_keys:
        if 'platform' not in results[key][fastq].keys():
          results[key][fastq]['platform'] = 'illumina'
          ill_fastq.append(fastq)
          results[key]['illumina'] = { 'platform' : 'illumina' }

      if len(ill_fastq) > 1 :
        for seqkit_val in ['num_seqs', 'sum_len', 'sum_gap', 'N50_num', 'sum_n']:
          results[key]['illumina'][seqkit_val] = round(float(results[key][ill_fastq[0]][seqkit_val]) + float(results[key][ill_fastq[1]][seqkit_val]), 0)

        for seqkit_val in ['avg_len', 'Q1', 'Q2', 'Q3', 'N50', 'Q20(%)', 'Q30(%)', 'AvgQual', 'GC(%)']:
          results[key]['illumina'][seqkit_val] = round((float(results[key][ill_fastq[0]][seqkit_val]) + float(results[key][ill_fastq[1]][seqkit_val]))/2, 2)

        results[key]['illumina']['file'] = f"{results[key][ill_fastq[0]]['file']},{results[key][ill_fastq[1]]['file']}"

        results[key]['illumina']['min_len'] = min(
          float(results[key][ill_fastq[0]]['min_len']), 
          float(results[key][ill_fastq[1]]['min_len']))
        
        results[key]['illumina']['max_len'] = max(
          float(results[key][ill_fastq[0]]['max_len']), 
          float(results[key][ill_fastq[1]]['max_len']))
        
        del results[key][ill_fastq[0]]
        del results[key][ill_fastq[1]]
    return results

def tsv_file(results_dict):
    final_results_dict = {}
    
    # converting the final dict to something tsv friendly
    sorted_keys = list(sorted(results_dict.keys()))
    all_keys = []
    print(sorted_keys)
    for sample in sorted_keys:
      final_results_dict[sample] = {}
      if 'nanopore' in results_dict[sample].keys():
        for result in results_dict[sample]['nanopore'].keys():
          final_results_dict[sample][f"seqkit_{result}"] = results_dict[sample]['nanopore'][result]

      if 'illumina' in results_dict[sample].keys():
        for result in results_dict[sample]['illumina'].keys():
          final_results_dict[sample][f"seqkit_illumina_{result}"] = results_dict[sample]['illumina'][result]

      if 'mash' in results_dict[sample].keys():
        for result in results_dict[sample]['mash'].keys():
          final_results_dict[sample][f"mash_{result}"] = results_dict[sample]['mash'][result]

      for analysis in ['pypolca', 'busco']:
        if analysis in results_dict[sample].keys():
          for assembler in results_dict[sample][analysis].keys():
            for result in results_dict[sample][analysis][assembler].keys():
              final_results_dict[sample][f"{assembler}_{analysis}_{result}"] = results_dict[sample][analysis][assembler][result]

      if 'assembly_info' in results_dict[sample].keys():
        for assembler in results_dict[sample]['assembly_info'].keys():
          for result in ['num_contigs', 'total_length', 'circ_contigs']:
            final_results_dict[sample][f"{assembler}_{result}"] = results_dict[sample]['assembly_info'][assembler][result]

      if 'circulocov' in results_dict[sample].keys():
        for assembler in results_dict[sample]['circulocov'].keys():
          for result in results_dict[sample]['circulocov'][assembler]['all'].keys():
            final_results_dict[sample][f"{assembler}_circulocov_{result}"] = results_dict[sample]['circulocov'][assembler]['all'][result]
      
      all_keys += list(final_results_dict[sample].keys())

    unique_key = list(set(all_keys))

    with open('donut_falls_summary.tsv', 'w') as tsv:
      i = 0
      for sample in sorted_keys:
        for key in unique_key:
          if key not in final_results_dict[sample].keys():
            final_results_dict[sample][key] = ""

        final_results_dict[sample] = { "sample": sample, **dict(sorted(final_results_dict[sample].items()))}

        possible_fieldnames = [
            "sample",
            "seqkit_num_seqs",
            "seqkit_avg_len",
            "seqkit_AvgQual",
            "seqkit_GC(%)",
            "mash_dist",
            "flye_total_length",
            "flye_num_contigs",
            "flye_circ_contigs",
            "flye_circulocov_nanopore_meandepth",
            "flye_circulocov_unmapped_nanopore_pc",
            "flye_circulocov_illumina_meandepth",
            "flye_circulocov_unmapped_illumina_pc",
            "flye_busco_reoriented",
            "flye_busco_clair3",
            "flye_busco_polypolish",
            "flye_busco_pypolca",
            "flye_pypolca_Insertion/Deletion_Errors_Found",
            "flye_pypolca_Substitution_Errors_Found",
            "myloasm_total_length",
            "myloasm_num_contigs",
            "myloasm_circ_contigs",
            "myloasm_circulocov_nanopore_meandepth",
            "myloasm_circulocov_unmapped_nanopore_pc",
            "myloasm_circulocov_illumina_meandepth",
            "myloasm_circulocov_unmapped_illumina_pc",
            "myloasm_busco_reoriented",
            "myloasm_busco_clair3",
            "myloasm_busco_polypolish",
            "myloasm_busco_pypolca",
            "myloasm_pypolca_Insertion/Deletion_Errors_Found",
            "myloasm_pypolca_Substitution_Errors_Found",
            "raven_total_length",
            "raven_num_contigs",
            "raven_circ_contigs",
            "raven_circulocov_nanopore_meandepth",
            "raven_circulocov_unmapped_nanopore_pc",
            "raven_circulocov_illumina_meandepth",
            "raven_circulocov_unmapped_illumina_pc",
            "raven_busco_reoriented",
            "raven_busco_clair3",
            "raven_busco_polypolish",
            "raven_busco_pypolca",
            "raven_pypolca_Insertion/Deletion_Errors_Found",
            "raven_pypolca_Substitution_Errors_Found",
            "unicycler_total_length",
            "unicycler_num_contigs",
            "unicycler_circ_contigs",
            "unicycler_circulocov_nanopore_meandepth",
            "unicycler_circulocov_unmapped_nanopore_pc",
            "unicycler_circulocov_illumina_meandepth",
            "unicycler_circulocov_unmapped_illumina_pc",
            "unicycler_busco_unicycler"
          ]

        fieldnames = []
        for name in possible_fieldnames:
          if name in final_results_dict[sample].keys():
            fieldnames.append(name)

        w = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
        if i < 1 :
          w.writeheader()
          i += 1
        filtered_row = {k: final_results_dict[sample][k] for k in fieldnames}
        w.writerow(filtered_row)

seqkit_dict     = seqkit_file('seqkit_summary.tsv') if exists('seqkit_summary.tsv') else {}
    
if not seqkit_dict:
  print('FATAL : Something is wrong and seqkit results were not located.')
  exit(1)

mash_dict          = mash_file('mash_summary.tsv') if exists('mash_summary.tsv') else {}
pypolca_dict       = pypolca_file('pypolca_summary.tsv') if exists('pypolca_summary.tsv') else {}
assembly_info_dict = assembly_info_file('assembly_info.csv') if exists('assembly_info.csv') else {}
busco_dict         = busco_results()
circulocov_dict    = circulocov_results()

final_results      = combine_results(seqkit_dict, mash_dict, pypolca_dict, assembly_info_dict, busco_dict, circulocov_dict)

final_file(final_results)
tsv_file(final_results)
