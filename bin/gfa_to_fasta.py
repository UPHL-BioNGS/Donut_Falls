#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

def convert_to_fasta(records, gfa_file):
    gfa_path = Path(gfa_file)
    outfile = gfa_path.with_suffix(".fasta")

    # Convert list of records to a dict for fast lookup
    rec_dict = {r['seq_name']: r for r in records}

    with open(gfa_file, 'r') as file, open(outfile, 'w') as output_file:
        for line in file:
            parts = line.strip().split()
            if parts and parts[0] == "S":
                header = parts[1]
                seq = parts[2]
                if header in rec_dict:
                    rec = rec_dict[header]
                    new_header = f">{header} length={rec['length']} circular={rec['circ.']} gc_per={rec['GC content %']}\n"
                    output_file.write(new_header)
                    output_file.write(seq + "\n")


def read_gfastats(gfastats_file):
    records = []  # must be a list to append dicts
    with open(gfastats_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("Seq"):  # skip header or empty lines
                continue
            parts = line.split()
            records.append({
                "sample": "prefix",
                "assembler": "raven",
                "seq_name": parts[1],
                "length": parts[3],
                "cov.": "",
                "circ.": parts[-1],
                "repeat": "",
                "mult.": "",
                "alt_group": "",
                "graph_path": "",
                "GC content %": parts[-3]  # assuming GC content is column 8 (0-based indexing)
            })
    return records

# Read stats and convert
records = read_gfastats("df_raven_gfastats.txt")
convert_to_fasta(records, "df_raven.gfa")

# Create DataFrame and save selected columns
df = pd.DataFrame(records)
print(df.to_csv(sep="\t", index=False))
df.to_csv(
    "prefix_assembly_info.txt",
    sep="\t",
    columns=["sample", "assembler", "seq_name", "length", "cov.", "circ.", "repeat", "mult.", "alt_group", "graph_path"],
    index=False
)
