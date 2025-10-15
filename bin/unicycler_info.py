#!/usr/bin/env python3
import pandas as pd

records = []
with open("2114703_unicycler.fasta") as f:
  for line in f:
    if line.startswith(">"):
      line = line.strip().replace(">", "")
      parts = line.split()
      circular="N"
      if "circular=true" in line:
        circular="Y"

      records.append({
        "sample": "prefix",
        "assembler": "unicycler",
        "seq_name": parts[0],
        "length": parts[1].split("=")[1],
        "cov.": parts[2].split("=")[1],
        "circ.": circular,
        "repeat": "",
        "mult.": "",
        "alt_group": "",
        "graph_path": ""
      })
    
df = pd.DataFrame(records)
df.to_csv("prefix_unicycler_assembly_info.txt", sep="\t", index=False)
