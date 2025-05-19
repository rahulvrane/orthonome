import re
import sys
from tqdm import tqdm

# Function to convert the original headers to the simplified format
def convert_headers(header):
    protein_id_match = re.search(r"\[protein_id=([XY]P_\d+\.\d+)\]", header)
    gene_match = re.search(r"\[gene=([^\]]+)\]", header)
    
    if protein_id_match and gene_match:
        protein_id = protein_id_match.group(1)
        gene = gene_match.group(1)
        return f">{protein_id}\t{gene}"
    else:
        return None

# Function to process a FASTA file
def process_fasta_file(input_file, output_file):
    total_genes = 0
    xp_genes = 0
    yp_genes = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in tqdm(infile):
            if line.startswith('>'):
                modified_header = convert_headers(line.strip())
#                print(modified_header)
                if modified_header:
                    outfile.write(modified_header + '\n')
                    total_genes += 1
                    
                    # Count genes by prefix
                    if "XP_" in modified_header:
                        xp_genes += 1
                    elif "YP_" in modified_header:
                        yp_genes += 1
            else:
                outfile.write(line)
    
    # Print gene counts
    print(f"\nGene counts:")
    print(f"Total genes processed: {total_genes}")
    print(f"Nuclear genes (XP_): {xp_genes}")
    print(f"Mitochondrial genes (YP_): {yp_genes}")

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py <input_fasta_file> <output_fasta_file>")
    sys.exit(1)

# Get the input and output file paths from command line arguments
input_fasta_file = sys.argv[1]
output_fasta_file = sys.argv[2]

# Process the input FASTA file and write to the output file
process_fasta_file(input_fasta_file, output_fasta_file)
