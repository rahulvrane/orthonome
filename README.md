# README #

Please cite: Rane, RV, Oakeshott, JG, Nguyen, T., Hoffmann, AA, & Lee, SF (2017). Ortho-a new pipeline for predicting high quality orthologous gene sets applicable to complete and draft genomes. BMC genomics , 18 (1), 673. 

https://doi.org/10.1186/s12864-017-4079-6

For searchable Drosophila orthologue predictions for 20 Flybase + modENCODE genomes please visit http://www.orthonome.com 


### Who do I talk to? ###

* orthonome@gmail.com

##########################################################################

# Orthonome Python 3 CLI Usage #

The Orthonome pipeline is now executed using the `orthonome_cli.py` Python 3 script. This command-line interface (CLI) provides a structured way to run the entire pipeline or individual steps.

**Basic Usage:**
```bash
python /path/to/orthonome_cli.py [OPTIONS] COMMAND [ARGS]...
```
Or, if `orthonome_cli.py` is in your current directory and executable:
```bash
./orthonome_cli.py [OPTIONS] COMMAND [ARGS]...
```

**Getting Help:**
To see the main help message, including all available subcommands:
```bash
python /path/to/orthonome_cli.py --help
```
To get help for a specific subcommand (you may need to provide required global options like `--species_prefix_file`):
```bash
python /path/to/orthonome_cli.py --species_prefix_file your_species_list.txt <subcommand> --help
```

**Global Options:**
These options apply to the main command and are used by various subcommands:
*   `--species_prefix_file FILE`: Path to the file containing species prefixes, one per line. (Required)
*   `--threads INTEGER`: Number of threads to use for parallelizable tasks. (Required, default: 4)
*   `--orthonome_dir DIRECTORY`: Path to the Orthonome installation directory. If not provided, the CLI attempts to auto-detect it (assuming `orthonome_cli.py` is in the Orthonome root directory).

**Available Subcommands:**
*   `prepare-inputs`: Runs `modify_info.py` to process `.preinfo` files into `.info` files and validates that all required input files (`.pep`, `.nuc`, `.info`) for each species prefix are present and not empty.
*   `run-diamond-blast`: Runs DIAMOND `makedb` for all species, then DIAMOND `blastp` for all species pairs, and finally `mcl` clustering on the blastp results.
*   `run-sw-align`: Runs Smith-Waterman alignments for all species pairs using `sw_runner.sh` and `ssearch36`.
*   `process-pairwise-comparisons`: Generates necessary intermediate files (`genelists.txt`, `Spp_list.idx`, `combinations.txt`) and then creates and runs `PairComparisonN.sh` scripts for each species pair to produce MSOAR inputs.
*   `generate-phylogeny`: Generates a consensus phylogeny using `Ortholog_pairs_to_FastTreephy.py` and prepares the tree file for MultiMSOAR.
*   `run-multimsoar-and-summarize`: Runs `MultiMSOAR` to identify orthogroups and then uses `summarise_orthogroups_internet_OUT.py` to produce final summary files.
*   `run-all`: Runs the complete Orthonome pipeline by executing all the above subcommands in the correct sequence.

**Example: Running the Full Pipeline**
```bash
python /path/to/orthonome/orthonome_cli.py --species_prefix_file /data/my_species_list.txt --threads 8 run-all
```

##########################################################################

# Installation #

To get started, run the following commands in the directory where you want to install Orthonome. 

## Clone from github ##
```
git clone https://github.com/rahulvrane/orthonome.git
cd orthonome
```
## OR Get the latest release ##
```
wget www.orthonome.com/help/orthonome_current_release.zip
unzip orthonome_current_release.zip
cd orthonome
```
## Make scripts and programs executable ##
It's recommended to ensure all scripts are executable:
```
chmod +x *.py *.sh
chmod +x Programs/*sh Programs/*py
# Compile C/C++ programs
cd Programs/
make
cd ..
# Compile phylip (if needed, though not directly used by CLI)
# cd Programs/phylip-3.67/
# make install
# cd ../..
```
Note: The Python scripts have been updated to Python 3 and are now orchestrated by `orthonome_cli.py`. The C/C++ programs in the `Programs/` directory are still used.

##########################################################################

# Dependencies #

Orthonome relies on several external command-line tools and Python libraries.

**Command-Line Tools:**
The `orthonome_cli.py` script will check for the following required tools in your system's PATH:
*   `mafft`: For multiple sequence alignment.
*   `gffread`: For GFF3 file processing.
*   `mcl`: For MCL clustering.
*   `diamond`: For protein sequence alignment.
*   `ssearch36`: (from the FASTA toolkit) For Smith-Waterman alignments.
*   `parallel`: (GNU Parallel) For parallelizing tasks.

Please ensure these are installed and accessible. The original Orthonome also listed FastTree2 and OpenMPI; these might be used by scripts called within the pipeline (e.g., `Ortholog_pairs_to_FastTreephy.py`) but are not directly checked by the main CLI script. NCBI BLAST+ is also a common dependency in bioinformatics, though DIAMOND is used as the primary aligner in the CLI workflow.

**Python Environment:**
*   **Python 3:** All scripts in Orthonome have been upgraded to Python 3 (Python 3.7+ recommended).
*   **Python Libraries:** Key Python libraries used by Orthonome and its CLI include:
    *   `click` (for the command-line interface)
    *   `numpy`
    *   `natsort`
    *   `tqdm`
    *   `biopython` (specifically `Bio.Seq` for translation tasks)

It is highly recommended to use a Python virtual environment (e.g., using `venv` or `conda`) to manage these dependencies. You can install them using pip:
```bash
pip install click numpy natsort tqdm biopython
```

##########################################################################

# Inputs files and preparing run inputs #

## Required inputs: ##
For every species to be analysed, the pipeline requires three input files starting with the same prefix and with the indicated suffixes

1. Gene nucleotide sequences (`<prefix>.nuc` or `<prefix>.cds` - the CLI expects `.nuc`)
2. Gene protein sequences (`<prefix>.pep`)
3. Gene information file (`<prefix>.info`). This is typically generated from a GFF3 file using the `prepare-inputs` step or the older `orthonome_inputs_after_gffread.py` script.

Additionally, the pipeline also needs a file listing all prefixes (e.g., `species_list.txt`, provided via `--species_prefix_file` to the CLI).

## If starting from a genome assembly and annotation in GFF3 format (NCBI format)

The genome and gff3 formats can be converted into the above inputs. This requires the GFFREAD utility. The `orthonome_inputs_after_gffread.py` script can assist with this, or you can use the `gffread` utility directly.

**Using `orthonome_inputs_after_gffread.py` (Python 3 compatible):**
```bash
# PRF should be your species prefix
# PRF_ncbi.gff is your input GFF3 file
# PRF_genome.fasta is your input genome FASTA file

# Extract all cds sequences
gffread PRF_ncbi.gff -g PRF_genome.fasta -x PRF.cds

# Extract one transcript per gene if multiple isoforms have been annotated
grep ">" PRF.cds | tr ' ' '\t' | sort -k2,2 -u | cut -f1 | tr -d ">" > PRF.idx

# Convert the genomic coordinates and cds sequences to the required formats
python /path/to/orthonome/orthonome_inputs_after_gffread.py PRF PRF_ncbi.gff
```
This will produce `PRF.pep`, `PRF.nuc` (from `PRF.cds`), and `PRF.preinfo`. The `PRF.preinfo` file must then be processed by the `prepare-inputs` step of `orthonome_cli.py` (which runs `modify_info.py`) to produce `PRF.info`.

**Alternative GFF processing (manual steps):**
In some cases - if the GFF has a non-NCBI-like format - the script `orthonome_inputs_after_gffread.py` might need adjustment or manual steps might be preferred:
```bash
CODE=PRF
GFF_FILE=your_annotation.gff3
GENOME_FILE=your_genome.fasta

gffread "$GFF_FILE" -g "$GENOME_FILE" -x "$CODE".cds -y "$CODE".prot 2>/dev/null && echo "GFF converted"
grep ">" "$CODE".prot | tr ' ' '\t' | sort -k2,2 -u | tr -d ">" > "$CODE".selected_genes.IDX
# Ensure seqtk is installed if using this part of the workflow
# seqtk subseq "$CODE".prot $CODE.selected_genes.IDX | sed 's/>.*=/>/g' > "$CODE".pep && echo "PEP created"
# seqtk subseq "$CODE".cds $CODE.selected_genes.IDX | sed 's/>.*=/>/g' > "$CODE".nuc && echo "NUC created"

# If not using seqtk, ensure .pep and .nuc files are generated by other means from selected transcripts.
# The .pep and .nuc files should correspond to the entries in $CODE.selected_genes.IDX.

awk '$3=="gene"' "$GFF_FILE" | sed 's/ID=.*Name=//g' | cut -f1 -d ';' | awk '{FS=OFS="\t"}{print $9,$9,$1,$7,$4}' > "$CODE".preinfo

# Then run the CLI:
# python /path/to/orthonome_cli.py --species_prefix_file species_list.txt prepare-inputs
# This will process $CODE.preinfo into $CODE.info.
```
Ensure the final `<prefix>.pep`, `<prefix>.nuc`, and `<prefix>.info` files are in your working directory before running subsequent Orthonome CLI steps.

##########################################################################

# Python 3 Upgrade Note #

All Python scripts included in the Orthonome suite have been upgraded from Python 2 to Python 3. The pipeline is now orchestrated via the new Python 3 based command-line interface `orthonome_cli.py`. The original `Orthonome_run.sh` script is no longer the primary method for running the pipeline, though its internal logic was used as a reference for the CLI subcommands.

##########################################################################
