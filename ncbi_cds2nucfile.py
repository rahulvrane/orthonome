import re
import sys
from tqdm import tqdm

# Function to convert the original headers to the simplified format
def convert_headers(header):
    protein_id_match = re.search(r"\[protein_id=(XP_\d+\.\d+)\]", header)
    gene_match = re.search(r"\[gene=([^\]]+)\]", header)
    
    if protein_id_match and gene_match:
        protein_id = protein_id_match.group(1)
        gene = gene_match.group(1)
        return f">{protein_id}\t{gene}"
    else:
        return None

# Function to process a FASTA file
def process_fasta_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in tqdm(infile):
            if line.startswith('>'):
                modified_header = convert_headers(line.strip())
#                print(modified_header)
                if modified_header:
                    outfile.write(modified_header + '\n')
            else:
                outfile.write(line)

# Check if the correct number of arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py <input_fasta_file> <output_fasta_file>")
    sys.exit(1)

# Get the input and output file paths from command line arguments
input_fasta_file = sys.argv[1]
output_fasta_file = sys.argv[2]

# Process the input FASTA file and write to the output file
process_fasta_file(input_fasta_file, output_fasta_file)
