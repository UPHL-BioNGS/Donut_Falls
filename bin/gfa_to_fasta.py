import csv
import glob

def convert_to_fasta(summary_dict, gfa_file):
    outfile = '_'.join(gfa_file.split('.')[:-1]) + ".fasta"
    with open(gfa_file, mode='r') as file:
        for line in file:
            parts = line.split()
            if parts and parts[0] == "S":
                header = parts[1]
                seq = parts[2]
                if header in summary_dict.keys():
                    new_header = ">" + header + " length=" + summary_dict[header]['Total segment length'] + " circular=" + summary_dict[header]["circular"].replace("N","false").replace("Y","true") + " gc_per=" + summary_dict[header]["GC content %"] + "\n"
                    with open(outfile, mode='a') as output_file:
                        output_file.write(new_header)
                        output_file.write(seq + "\n")

def read_summary_csv(gfastats_file):
    summary_dict = {}
    with open(gfastats_file, mode='r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            key = row['Header']
            summary_dict[key] = row
            with open("noncircular.txt", mode='a') as output_file:
                if summary_dict[key]["circular"] == "N":
                    output_file.write(key + "\n")
    return summary_dict

gfastats_file = glob.glob("*_gfastats_summary.csv")
gfa_file = glob.glob("*.gfa")

summary_dict = read_summary_csv(gfastats_file[0])
convert_to_fasta(summary_dict, gfa_file[0])
