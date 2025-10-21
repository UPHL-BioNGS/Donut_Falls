#!/usr/bin/env python3
import re
import pandas as pd

records = []
with open("assembly_primary.fa") as f:
    for line in f:
        if line.startswith(">"):
            line = line.strip().replace(">", "")
            parts = re.split(r"[ _]", line)

            records.append({
                "sample": "prefix",
                "assembler": "myloasm",
                "seq_name": parts[0],
                "length": parts[1].split("-")[1],
                "cov.": parts[3].split("-")[1],
                "circ.": {"yes": "Y", "no": "N", "possibly": "N"}.get(parts[2].split("-")[1], "N"),
                "repeat": {"yes": "Y", "no": "N", "possibly": "N"}.get(parts[4].split("-")[1], "N"),
                "mult.": parts[5].split("=")[1],
                "alt_group": "",
                "graph_path": ""
            })
    
df = pd.DataFrame(records)
df.to_csv("prefix_assembly_info.txt", sep="\t", index=False)

