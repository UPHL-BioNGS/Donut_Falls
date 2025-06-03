#!/usr/bin/env python3
import glob
import json
import csv
from os.path import exists


def busco_results():
    results = {}
    files = glob.glob('short_summary*txt')
    for file in files:
      sample, assembler, step = file.split('.')[-2].rsplit("_", 2)
      if sample not in results.keys():
        results[sample] = {}
      if assembler not in results[sample].keys():
        results[sample][assembler] = {}      
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
      sample, assembler = file.replace('_reoriented_overall_summary.txt','').rsplit("_", 1)
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

def combine_results(seqkit_dict, mash_dict, pypolca_dict, gfastats_dict, busco_dict, circulocov_dict):
      final_results = seqkit_dict

      for key in final_results.keys():

        if key in mash_dict.keys():
          final_results[key]['mash'] = mash_dict[key]

        if key in pypolca_dict.keys():
          final_results[key]['pypolca'] = pypolca_dict[key]

        if key in gfastats_dict.keys():
          final_results[key]['gfastats'] = gfastats_dict[key]

        if key in busco_dict.keys():
          final_results[key]['busco'] = busco_dict[key]

        if key in circulocov_dict.keys():
          final_results[key]['circulocov'] = circulocov_dict[key]

        #final_results[key]['assemblers'] = '${params.assembler}'
        final_results[key]['assemblers'] = 'flye,unicycler,raven,myloasm'
      
      return final_results
        
def final_file(results):
    with open('donut_falls_summary.json', 'w') as json_file:
      json.dump(results, json_file, indent=4)

def gfastats_file(file):
    results = {}
    with open(file, mode='r', newline='') as file:
      reader = csv.DictReader(file, delimiter=",")
      for row in reader:
        sample, assembler = row['sample'].rsplit("_", 1)
        if sample not in results.keys():
          results[sample] = {}
        if assembler not in results[sample].keys():
          results[sample][assembler] = {}
        results[sample][assembler][row['Header']] = row

    for sample in results.keys():
      for assembler in results[sample].keys():
        total_length = 0
        num_circular = 0
        contigs = list(results[sample][assembler].keys())
        results[sample][assembler]['num_contigs'] = len(contigs)
        for contig in contigs:
          total_length = total_length + int(results[sample][assembler][contig]['Total segment length'])
          if results[sample][assembler][contig]['Is circular'] == 'Y':
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
    with open('donut_falls_summary.tsv', 'w') as tsv:
      i = 0
      sorted_keys = list(sorted(results_dict.keys()))
      for sample in sorted_keys:
        final_results_dict[sample] = {}
        for result_key in ['sample', 'num_seqs', 'avg_len', 'AvgQual']:
          final_results_dict[sample][result_key] = results_dict[sample]['nanopore'][result_key] if results_dict[sample]['nanopore'][result_key] else 0
        
        final_results_dict[sample]['assemblers'] = results_dict[sample]['assemblers']
        final_results_dict[sample]['total_illumina_reads'] = results_dict[sample]['illumina']['num_seqs'] if results_dict[sample]['illumina']['num_seqs'] else 0
        final_results_dict[sample]['nanopore_illumina_mash_distance'] = results_dict[sample]['mash']['dist'] if  results_dict[sample]['mash']['dist'] else 1

        for assembler in ['flye', 'raven', 'myloasm']:
          if assembler in results_dict[sample]['assemblers']:
            if assembler in list(results_dict[sample]['gfastats'].keys()):
              for result in ['total_length', 'num_contigs', 'circ_contigs']:
                final_results_dict[sample][f"{assembler}_{result}"] = results_dict[sample]['gfastats'][assembler][result] if results_dict[sample]['gfastats'][assembler][result] else 0
              for result in ['nanopore_meandepth', 'nanopore_meandepth', 'unmapped_nanopore', 'unmapped_nanopore_pc', 'illumina_meandepth', 'unmapped_illumina', 'unmapped_illumina_pc', 'illumina_meandepth']:
                final_results_dict[sample][f"{assembler}_{result}"] = results_dict[sample]['circulocov'][assembler]['all'][result] if results_dict[sample]['circulocov'][assembler]['all'][result] else 0
              print(results_dict[sample]['pypolca'][assembler])
              final_results_dict[sample][f"{assembler}_busco_polished"] = results_dict[sample]['busco'][assembler]['pypolca']
              final_results_dict[sample][f"{assembler}_quality_before_polishing"] = results_dict[sample]['pypolca'][assembler]['Consensus_Quality_Before_Polishing']
              final_results_dict[sample][f"{assembler}_QV_before_polishing"] = results_dict[sample]['pypolca'][assembler]['Consensus_QV_Before_Polishing']
            else:
              for result in results + ['quality_before_polishing', 'QV_before_polishing' ]:
                final_results_dict[sample][assembler + '_' + result] = 0
              for result in ['busco', 'busco_polished']:
                final_results_dict[sample][assembler + '_' + result] = 'NF'
        if 'unicycler' in results_dict[sample]['assemblers']:
          if 'unicycler' in results_dict[sample].keys():
            for result in results + [ 'busco']:
              final_results_dict[sample]['unicycler_' + result] = results_dict[sample]['unicycler'][result]
          else:
            for result in results:
              final_results_dict[sample]['unicycler_' + result] = 0
            final_results_dict[sample]['unicycler_busco'] = 'NF'
        w = csv.DictWriter(tsv, final_results_dict[sample].keys(), delimiter='\t')
        if i < 1 :
          w.writeheader()
          i = i+1
        w.writerow(final_results_dict[sample])

def main():

    seqkit_dict     = seqkit_file('seqkit_summary.tsv') if exists('seqkit_summary.tsv') else {}
    
    if not seqkit_dict:
      print('FATAL : Something is wrong and seqkit results were not located.')
      exit(1)

    mash_dict       = mash_file('mash_summary.tsv') if exists('mash_summary.tsv') else {}
    pypolca_dict    = pypolca_file('pypolca_summary.tsv') if exists('pypolca_summary.tsv') else {}
    gfastats_dict   = gfastats_file('gfastats_summary.csv') if exists('gfastats_summary.csv') else {}
    busco_dict      = busco_results()
    circulocov_dict = circulocov_results()

    final_results   = combine_results(seqkit_dict, mash_dict, pypolca_dict, gfastats_dict, busco_dict, circulocov_dict)

    final_file(final_results)
    tsv_file(final_results)

if __name__ == '__main__':
  main()