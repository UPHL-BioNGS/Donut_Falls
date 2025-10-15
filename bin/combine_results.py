#!/usr/bin/env python3
"""
Donut Falls: Aggregate Sequencing QC and Assembly Summaries
Combines seqkit, mash, pypolca, assembly_info, BUSCO, and CirculoCOV results
into JSON and TSV summary outputs.
"""

import csv
import glob
import json
import logging
from pathlib import Path
from typing import Dict, Any

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)

# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def nested_dict_insert(results: Dict, keys: list[str], value: Any):
    """Insert a value into a nested dictionary given a list of keys."""
    d = results
    for key in keys[:-1]:
        d = d.setdefault(key, {})
    d[keys[-1]] = value

def safe_open_csv(file_path: str, delimiter: str = '\t') -> csv.DictReader | None:
    """Safely open a CSV/TSV file and return a DictReader."""
    try:
        return csv.DictReader(open(file_path, newline='', encoding='utf-8'), delimiter=delimiter)
    except FileNotFoundError:
        logging.warning(f"File not found: {file_path}")
        return None

# -----------------------------------------------------------------------------
# BUSCO
# -----------------------------------------------------------------------------
def busco_results() -> Dict[str, Any]:
    """Parse BUSCO short summary text files."""
    results: Dict[str, Any] = {}
    for file in glob.glob('short_summary*txt'):
        try:
            sample, assembler, step = Path(file).stem.split('.')[-1].rsplit("_", 2)
        except ValueError:
            logging.warning(f"Unexpected BUSCO filename format: {file}")
            continue
        with open(file, encoding='utf-8') as f:
            for line in f:
                if all(x in line for x in ['C:', 'S:', 'D:', 'F:', 'M:', 'n:']):
                    nested_dict_insert(results, [sample, assembler, step], line.strip())
                    break
    return results

# -----------------------------------------------------------------------------
# CirculoCOV
# -----------------------------------------------------------------------------
def circulocov_results() -> Dict[str, Any]:
    """Parse CirculoCOV overall summary files."""
    results: Dict[str, Any] = {}
    for file in glob.glob('*overall_summary.txt'):
        try:
            sample, assembler = Path(file).stem.replace('_reoriented_overall_summary', '').rsplit("_", 1)
        except ValueError:
            logging.warning(f"Unexpected CirculoCOV filename format: {file}")
            continue

        reader = safe_open_csv(file, delimiter='\t')
        if not reader:
            continue

        nested_dict_insert(results, [sample, assembler], {})

        for row in reader:
            results[sample][assembler][row['contigs']] = row

        all_stats = results[sample][assembler].setdefault('all', {'warnings': ""})
        missing = results[sample][assembler].get('missing', {})

        # Nanopore
        all_stats['unmapped_nanopore'] = missing.get('nanopore_numreads', 0)
        total_np = float(all_stats.get('nanopore_numreads', 1))
        all_stats['unmapped_nanopore_pc'] = round(float(all_stats['unmapped_nanopore']) / total_np, 2)
        if all_stats['unmapped_nanopore_pc'] > 0.1:
            all_stats['warnings'] += "High proportion of unmapped Nanopore reads,"

        # Illumina
        if 'illumina_numreads' in missing:
            all_stats['unmapped_illumina'] = missing['illumina_numreads']
            total_il = float(all_stats.get('illumina_numreads', 1))
            all_stats['unmapped_illumina_pc'] = round(float(all_stats['unmapped_illumina']) / total_il, 2)
            if all_stats['unmapped_illumina_pc'] > 0.1:
                all_stats['warnings'] += "High proportion of unmapped Illumina reads,"

    return results

# -----------------------------------------------------------------------------
# SeqKit
# -----------------------------------------------------------------------------
def seqkit_file(file_path: str) -> Dict[str, Any]:
    """Parse SeqKit summary TSV and merge Illumina pairs."""
    results: Dict[str, Any] = {}
    try:
        with open(file_path, newline='', encoding='utf-8') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            for row in reader:
                sample = row['sample']
                fastq_file = row['file']
                results.setdefault(sample, {})[fastq_file] = row
    except FileNotFoundError:
        logging.error(f"SeqKit file not found: {file_path}")
        return {}
    except KeyError as e:
        logging.error(f"SeqKit file missing expected column: {e}")
        return {}

    # Process each sample
    for sample, files in results.items():
        # Sort by avg_len descending: first entry is Nanopore
        sorted_files = dict(sorted(files.items(), key=lambda x: float(x[1].get('avg_len', 0)), reverse=True))
        results[sample] = sorted_files

        nanopore_file = next(iter(sorted_files))
        results[sample]['nanopore'] = sorted_files[nanopore_file]
        results[sample]['nanopore']['platform'] = 'nanopore'
        sorted_files.pop(nanopore_file)

        illumina_files = list(sorted_files.keys())
        results[sample]['illumina'] = {'platform': 'illumina'}

        if not illumina_files:
            logging.warning(f"No Illumina files detected for sample {sample}")
            continue

        if len(illumina_files) == 1:
            results[sample]['illumina'].update(sorted_files[illumina_files[0]])
        else:
            f1, f2 = illumina_files[:2]
            results[sample]['illumina']['file'] = f"{sorted_files[f1]['file']},{sorted_files[f2]['file']}"
            # Sum-based metrics
            sum_metrics = ['num_seqs', 'sum_len', 'sum_gap', 'N50_num', 'sum_n']
            for m in sum_metrics:
                val = float(sorted_files[f1].get(m, 0)) + float(sorted_files[f2].get(m, 0))
                results[sample]['illumina'][m] = round(val, 0)
            # Average-based metrics
            avg_metrics = ['avg_len', 'Q1', 'Q2', 'Q3', 'N50', 'Q20(%)', 'Q30(%)', 'AvgQual', 'GC(%)']
            for m in avg_metrics:
                v1 = float(sorted_files[f1].get(m, 0))
                v2 = float(sorted_files[f2].get(m, 0))
                results[sample]['illumina'][m] = round((v1 + v2)/2, 2)
            # Min/Max
            results[sample]['illumina']['min_len'] = min(float(sorted_files[f1].get('min_len', 0)), float(sorted_files[f2].get('min_len', 0)))
            results[sample]['illumina']['max_len'] = max(float(sorted_files[f1].get('max_len', 0)), float(sorted_files[f2].get('max_len', 0)))
        # Clean up
        for f in illumina_files:
            sorted_files.pop(f, None)

    logging.info(f"Parsed SeqKit results for {len(results)} samples from {file_path}")
    return results

# -----------------------------------------------------------------------------
# Mash
# -----------------------------------------------------------------------------
def mash_file(file_path: str) -> Dict[str, Any]:
    results: Dict[str, Any] = {}
    reader = safe_open_csv(file_path, delimiter='\t')
    if not reader:
        return results
    for row in reader:
        results[row['sample']] = row
    return results

# -----------------------------------------------------------------------------
# PyPolca
# -----------------------------------------------------------------------------
def pypolca_file(file_path: str) -> Dict[str, Any]:
    results: Dict[str, Any] = {}
    reader = safe_open_csv(file_path, delimiter='\t')
    if not reader:
        return results
    for row in reader:
        try:
            sample, assembler = row['sample'].rsplit("_", 1)
            results.setdefault(sample, {})[assembler] = row
        except ValueError:
            logging.warning(f"Skipping malformed sample in pypolca: {row['sample']}")
    return results

# -----------------------------------------------------------------------------
# Assembly Info
# -----------------------------------------------------------------------------
def assembly_info_file(file_path: str) -> Dict[str, Any]:
    results: Dict[str, Any] = {}
    reader = safe_open_csv(file_path, delimiter=',')
    if not reader:
        return results
    for row in reader:
        try:
            sample, assembler = row['sample'].rsplit("_", 1)
        except ValueError:
            logging.warning(f"Skipping malformed sample in assembly info: {row['sample']}")
            continue
        nested_dict_insert(results, [sample, assembler, row['Header']], row)

    for sample, assemblers in results.items():
        for assembler, contigs in assemblers.items():
            total_len = sum(int(c['Total segment length']) for c in contigs.values() if isinstance(c, dict))
            circ_count = sum(1 for c in contigs.values() if isinstance(c, dict) and c.get('Is circular') == 'Y')
            contigs['num_contigs'] = len(contigs)
            contigs['total_length'] = total_len
            contigs['circ_contigs'] = circ_count
    return results

# -----------------------------------------------------------------------------
# Combine Results
# -----------------------------------------------------------------------------
def combine_results(seqkit_dict, mash_dict, pypolca_dict, assembly_info_dict, busco_dict, circulocov_dict) -> Dict[str, Any]:
    final_results = seqkit_dict.copy()
    for sample in final_results:
        for name, source in [
            ('mash', mash_dict),
            ('pypolca', pypolca_dict),
            ('assembly_info', assembly_info_dict),
            ('busco', busco_dict),
            ('circulocov', circulocov_dict),
        ]:
            if sample in source:
                final_results[sample][name] = source[sample]
        final_results[sample]['assemblers'] = 'flye,unicycler,raven,myloasm'
    return final_results

# -----------------------------------------------------------------------------
# JSON Output
# -----------------------------------------------------------------------------
def write_json(results: Dict[str, Any], output_file: str = 'donut_falls_summary.json'):
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=4)
    logging.info(f"Wrote JSON summary: {output_file}")

# -----------------------------------------------------------------------------
# TSV Output
# -----------------------------------------------------------------------------
def tsv_file(results_dict: Dict[str, Any], output_file: str = 'donut_falls_summary.tsv'):
    logging.info(f"Writing flattened TSV summary: {output_file}")

    flattened_results: Dict[str, Any] = {}
    all_columns = set()

    def flatten_dict(d: Dict[str, Any], parent_key: str = "") -> Dict[str, Any]:
        items = {}
        for k, v in d.items():
            new_key = f"{parent_key}_{k}" if parent_key else k
            if isinstance(v, dict):
                items.update(flatten_dict(v, new_key))
            else:
                items[new_key] = v
        return items

    for sample, sample_data in results_dict.items():
        flat = flatten_dict(sample_data)
        flat['sample'] = sample
        flattened_results[sample] = flat
        all_columns.update(flat.keys())

    fieldnames = ['sample'] + sorted([c for c in all_columns if c != 'sample'])

    with open(output_file, 'w', newline='', encoding='utf-8') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for sample in sorted(flattened_results.keys()):
            row = {col: flattened_results[sample].get(col, "") for col in fieldnames}
            writer.writerow(row)

    logging.info(f"TSV summary written: {output_file}")

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    seqkit_path = 'seqkit_summary.tsv'
    if not Path(seqkit_path).exists():
        logging.error('FATAL: seqkit results not found.')
        exit(1)

    seqkit_dict = seqkit_file(seqkit_path)
    mash_dict = mash_file('mash_summary.tsv')
    pypolca_dict = pypolca_file('pypolca_summary.tsv')
    assembly_info_dict = assembly_info_file('assembly_info.csv')
    busco_dict = busco_results()
    circulocov_dict = circulocov_results()

    final_results = combine_results(seqkit_dict, mash_dict, pypolca_dict, assembly_info_dict, busco_dict, circulocov_dict)

    write_json(final_results)
    tsv_file(final_results)

    logging.info("Processing complete.")

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
