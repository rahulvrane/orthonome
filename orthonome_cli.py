#!/usr/bin/env python3
"""
Orthonome CLI - A command-line interface for Orthonome.
"""
import click
import os
import shutil
import sys
import glob
import subprocess
import itertools

REQUIRED_TOOLS = ["mafft", "gffread", "mcl", "diamond", "ssearch36", "parallel"]

def check_dependencies(dependencies: list[str]) -> None:
    """
    Checks if all required command-line tools are available in the system's PATH.
    Aborts with an error message if any dependency is missing.
    """
    for command in dependencies:
        if shutil.which(command) is None:
            click.echo(f"Error: Required command '{command}' not found in PATH.", file=sys.stderr)
            raise click.Abort()
    click.echo("All required dependencies found.")

@click.group()
@click.option('--species_prefix_file', 
              required=True, 
              type=click.Path(exists=True, dir_okay=False, readable=True),
              help="Path to the file containing species prefixes, one per line.")
@click.option('--threads', 
              required=True, 
              type=int, 
              default=4,
              show_default=True,
              help="Number of threads to use.")
@click.option('--orthonome_dir', 
              type=click.Path(file_okay=False, readable=True), 
              help="Path to the Orthonome installation directory. If not provided, attempts to auto-detect.")
def main_group(species_prefix_file, threads, orthonome_dir):
    """
    Orthonome CLI: A tool for running Orthonome workflows.
    """
    check_dependencies(REQUIRED_TOOLS)

    resolved_orthonome_dir = orthonome_dir
    if not resolved_orthonome_dir:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        resolved_orthonome_dir = script_dir 
        click.echo(f"Orthonome directory not provided, using auto-detected: {resolved_orthonome_dir}")
    
    if not os.path.isdir(resolved_orthonome_dir):
        click.echo(f"Error: Orthonome directory '{resolved_orthonome_dir}' does not exist or is not a directory.", file=sys.stderr)
        raise click.Abort()

    ctx = click.get_current_context()
    ctx.obj = {
        'species_prefix_file': species_prefix_file,
        'threads': threads,
        'orthonome_dir': resolved_orthonome_dir,
        'work_dir': os.getcwd()
    }
    click.echo(f"Main group initialized. Work dir: {ctx.obj['work_dir']}, Orthonome dir: {resolved_orthonome_dir}")

@main_group.command("prepare-inputs", help="Prepare and validate input files.")
@click.pass_context
def prepare_inputs_cmd(ctx):
    """
    Prepares input files by running modify_info.py and validates required files.
    """
    species_prefix_file = ctx.obj['species_prefix_file']
    orthonome_dir = ctx.obj['orthonome_dir']
    work_dir = ctx.obj['work_dir']

    click.echo("Starting input preparation...")

    modify_info_script_path = os.path.join(orthonome_dir, "modify_info.py")
    if not os.path.exists(modify_info_script_path):
        click.echo(f"Error: 'modify_info.py' not found at '{modify_info_script_path}'.", file=sys.stderr)
        raise click.Abort()
    if not os.access(modify_info_script_path, os.X_OK):
         click.echo(f"Warning: 'modify_info.py' at '{modify_info_script_path}' is not executable. Attempting to run with 'python'.")

    preinfo_files = glob.glob(os.path.join(work_dir, '*.preinfo'))
    if not preinfo_files:
        click.echo("No '*.preinfo' files found in the working directory. Skipping modify_info.py execution.")
    else:
        click.echo(f"Found {len(preinfo_files)} '.preinfo' files. Running modify_info.py...")
        for preinfo_file_path in preinfo_files:
            base_name = os.path.splitext(os.path.basename(preinfo_file_path))[0]
            info_file_path = os.path.join(work_dir, f"{base_name}.info")
            
            command = [sys.executable, modify_info_script_path, preinfo_file_path]
            
            click.echo(f"  Processing {preinfo_file_path} -> {info_file_path}")
            try:
                result = subprocess.run(command, capture_output=True, text=True, check=False)
                
                if result.returncode != 0:
                    click.echo(f"Error running modify_info.py for {preinfo_file_path}:", file=sys.stderr)
                    click.echo(f"  Stdout: {result.stdout}", file=sys.stderr)
                    click.echo(f"  Stderr: {result.stderr}", file=sys.stderr)
                    raise click.Abort()
                else:
                    with open(info_file_path, 'w') as outfile_succ:
                        outfile_succ.write(result.stdout)
            except Exception as e:
                click.echo(f"An exception occurred while running modify_info.py for {preinfo_file_path}: {e}", file=sys.stderr)
                raise click.Abort()
        click.echo("modify_info.py execution completed.")

    click.echo("Validating input files...")
    try:
        with open(species_prefix_file, 'r') as f:
            species_prefixes_list = [line.strip() for line in f if line.strip()] # Renamed to avoid conflict
    except IOError as e:
        click.echo(f"Error reading species prefix file '{species_prefix_file}': {e}", file=sys.stderr)
        raise click.Abort()

    if not species_prefixes_list:
        click.echo(f"Error: No species prefixes found in '{species_prefix_file}'.", file=sys.stderr)
        raise click.Abort()

    all_files_valid = True
    for prefix in species_prefixes_list:
        required_extensions = ['.pep', '.nuc', '.info']
        for ext in required_extensions:
            file_path = os.path.join(work_dir, f"{prefix}{ext}")
            if not os.path.exists(file_path):
                click.echo(f"Error: Required file '{file_path}' is missing for prefix '{prefix}'.", file=sys.stderr)
                all_files_valid = False
            elif os.path.getsize(file_path) == 0:
                click.echo(f"Error: Required file '{file_path}' is empty for prefix '{prefix}'.", file=sys.stderr)
                all_files_valid = False
    
    if not all_files_valid:
        raise click.Abort()

    click.echo("Input file validation successful. All required files are present and not empty.")
    click.echo("Input preparation and validation completed successfully.")

@main_group.command("run-diamond-blast", help="Run DIAMOND makedb, DIAMOND blastp, and MCL clustering.")
@click.pass_context
def run_diamond_blast_cmd(ctx):
    species_prefix_file_abs_path = os.path.abspath(ctx.obj['species_prefix_file'])
    threads = ctx.obj['threads']
    work_dir = ctx.obj['work_dir']
    
    blast_pairs_dir = os.path.join(work_dir, "blast_pairs")
    os.makedirs(blast_pairs_dir, exist_ok=True)
    click.echo(f"Created directory: {blast_pairs_dir}")

    try:
        with open(species_prefix_file_abs_path, 'r') as f:
            species_prefixes = [line.strip() for line in f if line.strip()]
    except IOError as e:
        click.echo(f"Error reading species prefix file '{species_prefix_file_abs_path}': {e}", file=sys.stderr)
        raise click.Abort()

    if not species_prefixes:
        click.echo(f"Error: No species prefixes found in '{species_prefix_file_abs_path}'.", file=sys.stderr)
        raise click.Abort()

    click.echo("Running DIAMOND makedb for all species...")
    makedb_command = f"parallel --gnu -j{threads} 'diamond makedb --in ../{{}}.pep --db {{}}' :::: {species_prefix_file_abs_path}"
    click.echo(f"Executing: {makedb_command} (cwd: {blast_pairs_dir})")
    try:
        result_makedb = subprocess.run(makedb_command, shell=True, check=True, cwd=blast_pairs_dir, 
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                       executable='/bin/bash') 
        click.echo("DIAMOND makedb completed successfully.")
        if result_makedb.stdout: click.echo(f"  Stdout:\n{result_makedb.stdout}")
        if result_makedb.stderr: click.echo(f"  Stderr:\n{result_makedb.stderr}")
    except subprocess.CalledProcessError as e:
        click.echo(f"Error during DIAMOND makedb: {e}", file=sys.stderr)
        raise click.Abort()

    click.echo("Running DIAMOND blastp for all species pairs...")
    species_pairs = list(itertools.permutations(species_prefixes, 2))
    
    for sp1, sp2 in species_pairs:
        success_file = os.path.join(blast_pairs_dir, f"{sp1}_{sp2}.success")
        if os.path.exists(success_file):
            click.echo(f"  Skipping {sp1} vs {sp2}, success file found: {success_file}")
            continue

        blastp_command = (
            f"parallel --gnu -j1 --bar " 
            f"'diamond blastp --query ../{{1}}.pep --db {{2}} --outfmt 6 --out {{1}}_{{2}}.blastp --threads {threads} && touch {{1}}_{{2}}.success' "
            f"::: {sp1} ::: {sp2}"
        )
        click.echo(f"  Running blastp for {sp1} vs {sp2}. Executing: {blastp_command} (cwd: {blast_pairs_dir})")
        try:
            result_blastp = subprocess.run(blastp_command, shell=True, check=True, cwd=blast_pairs_dir,
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                           executable='/bin/bash')
            if result_blastp.stdout: click.echo(f"    Stdout:\n{result_blastp.stdout}")
            if result_blastp.stderr: click.echo(f"    Stderr:\n{result_blastp.stderr}")
        except subprocess.CalledProcessError as e:
            click.echo(f"Error during DIAMOND blastp for {sp1} vs {sp2}: {e}", file=sys.stderr)
            raise click.Abort()

    click.echo("DIAMOND blastp for all pairs completed.")

    click.echo("Verifying all blastp jobs completed successfully...")
    all_blastp_successful = True
    for sp1, sp2 in species_pairs:
        success_file = os.path.join(blast_pairs_dir, f"{sp1}_{sp2}.success")
        if not os.path.exists(success_file):
            click.echo(f"Error: Blastp for {sp1} vs {sp2} did not complete successfully (missing {success_file}).", file=sys.stderr)
            all_blastp_successful = False
    
    if not all_blastp_successful:
        click.echo("One or more blastp jobs failed. Please check logs.", file=sys.stderr)
        raise click.Abort()
    click.echo("All blastp jobs verified.")

    click.echo("Running MCL clustering...")
    mcl_command = (
        f"cat *.blastp | grep -v \"#\" | "
        f"parallel --gnu --pipe -q awk '{{OFS=\"\\t\"}}{{if ($11<=0.5) print $1, $2, $12; else print $1,$2,0}}' | "
        f"mcl - --abc -q x -V all -te {threads} -o Allruns.clusters"
    )
    click.echo(f"Executing: {mcl_command} (cwd: {blast_pairs_dir})")
    try:
        result_mcl = subprocess.run(mcl_command, shell=True, check=True, cwd=blast_pairs_dir,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                    executable='/bin/bash')
        click.echo("MCL clustering completed successfully.")
        if result_mcl.stdout: click.echo(f"  Stdout:\n{result_mcl.stdout}")
        if result_mcl.stderr: click.echo(f"  Stderr:\n{result_mcl.stderr}")
    except subprocess.CalledProcessError as e:
        click.echo(f"Error during MCL clustering: {e}", file=sys.stderr)
        raise click.Abort()
    click.echo("run-diamond-blast subcommand completed successfully.")

@main_group.command("run-sw-align", help="Run Smith-Waterman alignments using sw_runner.sh.")
@click.pass_context
def run_sw_align_cmd(ctx):
    species_prefix_file_abs_path = os.path.abspath(ctx.obj['species_prefix_file'])
    threads = ctx.obj['threads'] 
    work_dir = ctx.obj['work_dir']
    orthonome_dir = ctx.obj['orthonome_dir']

    sw_scores_dir = os.path.join(work_dir, "sw_scores")
    os.makedirs(sw_scores_dir, exist_ok=True)
    sw_scores_dir_abs_path = os.path.abspath(sw_scores_dir)
    click.echo(f"Created directory: {sw_scores_dir_abs_path}")

    sw_runner_script_path = os.path.join(orthonome_dir, "sw_runner.sh")
    if not os.path.exists(sw_runner_script_path):
        click.echo(f"Error: 'sw_runner.sh' not found at '{sw_runner_script_path}'.", file=sys.stderr)
        raise click.Abort()
    if not os.access(sw_runner_script_path, os.X_OK):
        click.echo(f"Error: 'sw_runner.sh' at '{sw_runner_script_path}' is not executable.", file=sys.stderr)
        raise click.Abort()
    
    try:
        with open(species_prefix_file_abs_path, 'r') as f:
            species_prefixes = [line.strip() for line in f if line.strip()]
    except IOError as e:
        click.echo(f"Error reading species prefix file '{species_prefix_file_abs_path}': {e}", file=sys.stderr)
        raise click.Abort()

    if not species_prefixes:
        click.echo(f"Error: No species prefixes found in '{species_prefix_file_abs_path}'.", file=sys.stderr)
        raise click.Abort()

    click.echo("Running Smith-Waterman alignments for all species pairs...")
    species_pairs = list(itertools.permutations(species_prefixes, 2))
    
    commands_to_run = []
    for sp1, sp2 in species_pairs:
        success_file = os.path.join(sw_scores_dir_abs_path, f"{sp1}_{sp2}.success")
        if os.path.exists(success_file):
            click.echo(f"  Skipping {sp1} vs {sp2}, success file found: {success_file}")
            continue
        commands_to_run.append(f"bash {sw_runner_script_path} {sp1} {sp2} {threads} {sw_scores_dir_abs_path}")

    if commands_to_run:
        parallel_command_str = ( 
            f"parallel --gnu -j{threads} --bar --halt now,fail=1 " 
            "'{cmd}' ::: " + " ::: ".join([f'"{cmd}"' for cmd in commands_to_run])
        )
        if not commands_to_run: 
             parallel_command_str = "true" 

        click.echo(f"Executing {len(commands_to_run)} SW alignment jobs in parallel (up to {threads} concurrently)...")
        try:
            result_sw = subprocess.run(parallel_command_str, shell=True, check=True, cwd=work_dir,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                       executable='/bin/bash')
            if result_sw.stdout: click.echo(f"  Stdout (from parallel):\n{result_sw.stdout}")
            if result_sw.stderr: click.echo(f"  Stderr (from parallel):\n{result_sw.stderr}")
        except subprocess.CalledProcessError as e:
            click.echo(f"Error during Smith-Waterman alignments execution with parallel: {e}", file=sys.stderr)
            raise click.Abort()
    else:
        click.echo("All Smith-Waterman alignment pairs already processed (success files found).")

    click.echo("Smith-Waterman alignments for all pairs completed.")

    click.echo("Verifying all Smith-Waterman alignment jobs completed successfully...")
    all_sw_successful = True
    for sp1, sp2 in species_pairs:
        success_file = os.path.join(sw_scores_dir_abs_path, f"{sp1}_{sp2}.success")
        if not os.path.exists(success_file):
            click.echo(f"Error: Smith-Waterman alignment for {sp1} vs {sp2} did not complete successfully (missing {success_file}).", file=sys.stderr)
            all_sw_successful = False
    
    if not all_sw_successful:
        click.echo("One or more Smith-Waterman alignment jobs failed. Please check logs.", file=sys.stderr)
        raise click.Abort()
    click.echo("All Smith-Waterman alignment jobs verified.")
    click.echo("run-sw-align subcommand completed successfully.")

@main_group.command("process-pairwise-comparisons", help="Generate and run batch scripts for pairwise species comparisons.")
@click.pass_context
def process_pairwise_comparisons_cmd(ctx):
    species_prefix_file_abs_path = os.path.abspath(ctx.obj['species_prefix_file'])
    threads = ctx.obj['threads']
    work_dir = ctx.obj['work_dir']
    orthonome_dir = ctx.obj['orthonome_dir']

    click.echo("Starting processing of pairwise comparisons...")

    try:
        with open(species_prefix_file_abs_path, 'r') as f:
            species_prefixes = [line.strip() for line in f if line.strip()]
    except IOError as e:
        click.echo(f"Error reading species prefix file '{species_prefix_file_abs_path}': {e}", file=sys.stderr)
        raise click.Abort()
    if not species_prefixes:
        click.echo(f"Error: No species prefixes found in '{species_prefix_file_abs_path}'.", file=sys.stderr)
        raise click.Abort()

    genelists_file_path = os.path.join(work_dir, "genelists.txt")
    if os.path.exists(genelists_file_path): os.remove(genelists_file_path) 
    click.echo(f"Creating {genelists_file_path}...")
    for prefix in species_prefixes:
        info_file = os.path.join(work_dir, f"{prefix}.info")
        if not os.path.exists(info_file):
            click.echo(f"Error: Info file {info_file} not found for prefix {prefix}.", file=sys.stderr)
            raise click.Abort()
        try:
            with open(info_file, 'r') as f_info:
                gene_ids = [line.split('\t')[0] for line in f_info if line.strip()]
            with open(genelists_file_path, 'a') as f_genelists:
                f_genelists.write(f"{prefix}\t" + "\t".join(gene_ids) + "\n")
        except Exception as e:
            click.echo(f"Error processing info file {info_file} for genelists.txt: {e}", file=sys.stderr)
            raise click.Abort()
    click.echo("genelists.txt created successfully.")

    spp_list_idx_path = os.path.join(work_dir, "Spp_list.idx")
    click.echo(f"Creating {spp_list_idx_path}...")
    awk_command = f"awk '{{print $0\"_\"NR}}' \"{species_prefix_file_abs_path}\" > \"{spp_list_idx_path}\""
    try:
        subprocess.run(awk_command, shell=True, check=True, cwd=work_dir, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        click.echo(f"Error creating Spp_list.idx: {e}", file=sys.stderr)
        raise click.Abort()
    click.echo("Spp_list.idx created successfully.")

    combinations_script_path = os.path.join(orthonome_dir, "combinations.py")
    combinations_txt_path = os.path.join(work_dir, "combinations.txt")
    if not os.path.exists(combinations_script_path):
        click.echo(f"Error: 'combinations.py' not found at '{combinations_script_path}'.", file=sys.stderr)
        raise click.Abort()
    click.echo(f"Running combinations.py to create {combinations_txt_path}...")
    combinations_command = f"{sys.executable} \"{combinations_script_path}\" \"{spp_list_idx_path}\" | tr -d '() ' > \"{combinations_txt_path}\""
    try:
        subprocess.run(combinations_command, shell=True, check=True, cwd=work_dir, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        click.echo(f"Error running combinations.py: {e}", file=sys.stderr)
        raise click.Abort()
    click.echo("combinations.txt created successfully.")

    click.echo("Generating PairComparisonN.sh scripts...")
    pair_comparison_scripts = []
    combination_prf_entries = []
    process_species_pairs_script_path = os.path.join(orthonome_dir, "Process_Species_pairs.sh")
    if not os.path.exists(process_species_pairs_script_path):
        click.echo(f"Error: 'Process_Species_pairs.sh' not found at '{process_species_pairs_script_path}'.", file=sys.stderr)
        raise click.Abort()

    species_map = {str(i+1): prefix for i, prefix in enumerate(species_prefixes)}

    with open(combinations_txt_path, 'r') as f_comb:
        for i, line in enumerate(f_comb):
            p_idx_str, s_idx_str = line.strip().split(',')
            spp1 = species_map.get(p_idx_str)
            spp2 = species_map.get(s_idx_str)

            if not spp1 or not spp2:
                click.echo(f"Error: Invalid species indices {p_idx_str}, {s_idx_str} in combinations.txt", file=sys.stderr)
                raise click.Abort()

            pair_dir_name = f"{spp1}_{spp2}"
            combination_prf_entries.append(pair_dir_name)
            
            script_num = i + 1
            pair_comparison_script_path = os.path.join(work_dir, f"PairComparison{script_num}.sh")
            pair_comparison_scripts.append(pair_comparison_script_path)

            script_content = f"""#!/bin/bash
mkdir -p "{pair_dir_name}"
cd "{pair_dir_name}"
bash "{process_species_pairs_script_path}" \\
  "{os.path.join(work_dir, spp1 + '.pep')}" "{os.path.join(work_dir, spp2 + '.pep')}" \\
  "{os.path.join(work_dir, spp1 + '.nuc')}" "{os.path.join(work_dir, spp2 + '.nuc')}" \\
  "{os.path.join(work_dir, spp1 + '.info')}" "{os.path.join(work_dir, spp2 + '.info')}" \\
  "{spp1},{spp2}" "{p_idx_str}" "{s_idx_str}" 3 all \\
  1> PairComparison.log 2>&1 && touch "{os.path.join(work_dir, f'PairComparison{script_num}.success')}"
"""
            with open(pair_comparison_script_path, 'w') as f_script:
                f_script.write(script_content)
            os.chmod(pair_comparison_script_path, 0o755)
    
    combination_prf_path = os.path.join(work_dir, "combination_prf.txt")
    with open(combination_prf_path, 'w') as f_prf:
        for entry in combination_prf_entries:
            f_prf.write(entry + "\n")
    click.echo(f"Generated {len(pair_comparison_scripts)} PairComparisonN.sh scripts and combination_prf.txt.")

    msoar_inputs_dir = os.path.join(work_dir, "MultiMSOAR_inputs")
    os.makedirs(msoar_inputs_dir, exist_ok=True)
    click.echo(f"Created directory: {msoar_inputs_dir}")

    if pair_comparison_scripts:
        num_parallel_jobs = max(1, threads // 3) 
        click.echo(f"Running {len(pair_comparison_scripts)} PairComparison scripts in parallel (up to {num_parallel_jobs} concurrently)...")
        
        parallel_scripts_paths_str = ' '.join([f'"{s}"' for s in pair_comparison_scripts]) 
        parallel_command_str = f"parallel --gnu --bar -j{num_parallel_jobs} 'bash {{}}' ::: {parallel_scripts_paths_str}"
        
        click.echo(f"Executing: {parallel_command_str}")
        try:
            subprocess.run(parallel_command_str, shell=True, check=True, cwd=work_dir, executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            click.echo(f"Error running PairComparison scripts with parallel: {e}", file=sys.stderr)
            raise click.Abort()
        click.echo("PairComparison script execution completed.")
    else:
        click.echo("No PairComparison scripts to run.")

    click.echo("Verifying pairwise comparison results...")
    all_pairs_successful = True
    for i, script_path_iter in enumerate(pair_comparison_scripts): 
        script_num = i + 1
        success_marker = os.path.join(work_dir, f"PairComparison{script_num}.success")
        if not os.path.exists(success_marker):
            click.echo(f"Error: PairComparison{script_num}.sh (for {os.path.basename(script_path_iter)}) did not complete successfully (missing {success_marker}).", file=sys.stderr)
            all_pairs_successful = False
    
    with open(combination_prf_path, 'r') as f_prf:
        for prcomb in f_prf:
            prcomb = prcomb.strip()
            msoar_result_path = os.path.join(work_dir, prcomb, "MSOAR2_result")
            if not os.path.exists(msoar_result_path) or os.path.getsize(msoar_result_path) == 0:
                click.echo(f"Error: MSOAR2_result file is missing or empty for pair {prcomb} at {msoar_result_path}.", file=sys.stderr)
                all_pairs_successful = False
                
    if not all_pairs_successful:
        click.echo("One or more pairwise comparison jobs failed. Please check logs.", file=sys.stderr)
        raise click.Abort()
    
    click.echo("All pairwise comparison jobs verified successfully.")
    click.echo("process-pairwise-comparisons subcommand completed successfully.")

@main_group.command("generate-phylogeny", help="Generate consensus phylogeny and prepare tree for MultiMSOAR.")
@click.pass_context
def generate_phylogeny_cmd(ctx):
    species_prefix_file_abs_path = os.path.abspath(ctx.obj['species_prefix_file'])
    threads = ctx.obj['threads']
    work_dir = ctx.obj['work_dir']
    orthonome_dir = ctx.obj['orthonome_dir']

    click.echo("Starting phylogeny generation...")

    multi_msoar_inputs_dir = os.path.join(work_dir, "MultiMSOAR_inputs")
    if not os.path.isdir(multi_msoar_inputs_dir):
        click.echo(f"Error: Directory '{multi_msoar_inputs_dir}' not found. Please run 'process-pairwise-comparisons' first.", file=sys.stderr)
        raise click.Abort()

    script_path = os.path.join(orthonome_dir, "Ortholog_pairs_to_FastTreephy.py")
    if not os.path.exists(script_path):
        click.echo(f"Error: Script 'Ortholog_pairs_to_FastTreephy.py' not found at '{script_path}'.", file=sys.stderr)
        raise click.Abort()

    genelists_path = os.path.join(work_dir, "genelists.txt")
    
    try:
        with open(species_prefix_file_abs_path, 'r') as f:
            species_prefixes_for_outgroup = [line.strip() for line in f if line.strip()]
        if not species_prefixes_for_outgroup:
            click.echo(f"Error: No species found in '{species_prefix_file_abs_path}' to determine outgroup.", file=sys.stderr)
            raise click.Abort()
        outgroup_species = species_prefixes_for_outgroup[0] 
    except IOError as e:
        click.echo(f"Error reading species prefix file '{species_prefix_file_abs_path}': {e}", file=sys.stderr)
        raise click.Abort()

    command_phylogeny = [ # Renamed to avoid conflict
        sys.executable, script_path,
        "-l", genelists_path,
        "-g", multi_msoar_inputs_dir + "/", 
        "-p", work_dir, 
        "-n", work_dir, 
        "-t", str(threads),
        "--outgroup", outgroup_species 
    ]
    click.echo(f"Running Ortholog_pairs_to_FastTreephy.py in {multi_msoar_inputs_dir}...")
    click.echo(f"Command: {' '.join(command_phylogeny)}")
    try:
        result_fasttreephy = subprocess.run(command_phylogeny, cwd=multi_msoar_inputs_dir, check=True, 
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        click.echo("Ortholog_pairs_to_FastTreephy.py completed successfully.")
        if result_fasttreephy.stdout: click.echo(f"  Stdout:\n{result_fasttreephy.stdout}")
        if result_fasttreephy.stderr: click.echo(f"  Stderr:\n{result_fasttreephy.stderr}")
    except subprocess.CalledProcessError as e:
        click.echo(f"Error running Ortholog_pairs_to_FastTreephy.py: {e}", file=sys.stderr)
        click.echo(f"  Stdout: {e.stdout}", file=sys.stderr)
        click.echo(f"  Stderr: {e.stderr}", file=sys.stderr)
        raise click.Abort()

    nwk_file_path = os.path.join(multi_msoar_inputs_dir, "CONCAT_align_nuc.nwk")
    tree_dest_path = os.path.join(multi_msoar_inputs_dir, "Tree")

    if not os.path.exists(nwk_file_path) or os.path.getsize(nwk_file_path) == 0:
        click.echo(f"Error: Expected output file '{nwk_file_path}' not found or is empty.", file=sys.stderr)
        raise click.Abort()
    
    try:
        shutil.copy(nwk_file_path, tree_dest_path)
        click.echo(f"Copied '{nwk_file_path}' to '{tree_dest_path}'.")
    except IOError as e:
        click.echo(f"Error copying tree file: {e}", file=sys.stderr)
        raise click.Abort()

    click.echo(f"Modifying '{tree_dest_path}' file...")
    modify_tree_command = f"cat \"{species_prefix_file_abs_path}\" | awk '{{print $1 \"S\"(NR-1)}}' | while read line; do S=$(echo $line|cut -f1 -d ' '); PS=$(echo $line|cut -f2 -d ' '); sed -i \"s:$S:$PS:g\" Tree; done"
    click.echo(f"Executing: {modify_tree_command} (cwd: {multi_msoar_inputs_dir})")
    try:
        subprocess.run(modify_tree_command, shell=True, check=True, cwd=multi_msoar_inputs_dir, executable='/bin/bash')
        click.echo(f"Successfully modified '{tree_dest_path}'.")
    except subprocess.CalledProcessError as e:
        click.echo(f"Error modifying tree file '{tree_dest_path}': {e}", file=sys.stderr)
        raise click.Abort()
        
    click.echo("generate-phylogeny subcommand completed successfully.")

@main_group.command("run-multimsoar-and-summarize", help="Run MultiMSOAR and summarize orthogroups.")
@click.pass_context
def run_multimsoar_summarize_cmd(ctx):
    species_prefix_file_abs_path = os.path.abspath(ctx.obj['species_prefix_file'])
    work_dir = ctx.obj['work_dir']
    orthonome_dir = ctx.obj['orthonome_dir']
    # threads = ctx.obj['threads'] # Not directly used by MultiMSOAR or summarise script, but good to have if needed

    click.echo("Starting MultiMSOAR and summarization...")

    multi_msoar_inputs_dir = os.path.join(work_dir, "MultiMSOAR_inputs")
    blast_pairs_dir = os.path.join(work_dir, "blast_pairs") # For Allruns.clusters

    if not os.path.isdir(multi_msoar_inputs_dir):
        click.echo(f"Error: Directory '{multi_msoar_inputs_dir}' not found. Please run prerequisite steps.", file=sys.stderr)
        raise click.Abort()
    if not os.path.isdir(blast_pairs_dir):
        click.echo(f"Error: Directory '{blast_pairs_dir}' not found. Please run 'run-diamond-blast' first.", file=sys.stderr)
        raise click.Abort()

    # Run MultiMSOAR
    multimsoar_exe_path = os.path.join(orthonome_dir, "Programs", "MultiMSOAR") # Corrected path
    if not os.path.exists(multimsoar_exe_path):
        click.echo(f"Error: MultiMSOAR executable not found at '{multimsoar_exe_path}'.", file=sys.stderr)
        raise click.Abort()
    if not os.access(multimsoar_exe_path, os.X_OK):
        click.echo(f"Error: MultiMSOAR at '{multimsoar_exe_path}' is not executable.", file=sys.stderr)
        raise click.Abort()

    try:
        with open(species_prefix_file_abs_path, 'r') as f:
            num_species = len([line.strip() for line in f if line.strip()])
    except IOError as e:
        click.echo(f"Error reading species prefix file '{species_prefix_file_abs_path}': {e}", file=sys.stderr)
        raise click.Abort()
    
    if num_species == 0:
        click.echo(f"Error: No species found in '{species_prefix_file_abs_path}'. Cannot run MultiMSOAR.", file=sys.stderr)
        raise click.Abort()

    tree_file_path = os.path.join(multi_msoar_inputs_dir, "Tree")
    clusters_file_path = os.path.join(blast_pairs_dir, "Allruns.clusters")
    geneinfo_output_path = os.path.join(multi_msoar_inputs_dir, "Geneinfo") # Output file name
    orthogroups_output_path = os.path.join(multi_msoar_inputs_dir, "Orthogroups") # Output file name
    
    # Ensure prerequisite files for MultiMSOAR exist
    if not os.path.exists(tree_file_path):
        click.echo(f"Error: Tree file '{tree_file_path}' not found. Please run 'generate-phylogeny' first.", file=sys.stderr)
        raise click.Abort()
    if not os.path.exists(clusters_file_path):
        click.echo(f"Error: Clusters file '{clusters_file_path}' not found. Please run 'run-diamond-blast' first.", file=sys.stderr)
        raise click.Abort()


    multimsoar_command = [
        multimsoar_exe_path,
        str(num_species),
        "Tree", # Relative to CWD (multi_msoar_inputs_dir)
        os.path.relpath(clusters_file_path, multi_msoar_inputs_dir), # Relative path for clusters
        "Geneinfo", # Output, relative to CWD
        "Orthogroups" # Output, relative to CWD
    ]
    click.echo(f"Running MultiMSOAR in {multi_msoar_inputs_dir}...")
    click.echo(f"Command: {' '.join(multimsoar_command)}")
    try:
        result_msoar = subprocess.run(multimsoar_command, cwd=multi_msoar_inputs_dir, check=True,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        click.echo("MultiMSOAR completed successfully.")
        if result_msoar.stdout: click.echo(f"  Stdout:\n{result_msoar.stdout}")
        if result_msoar.stderr: click.echo(f"  Stderr:\n{result_msoar.stderr}")
    except subprocess.CalledProcessError as e:
        click.echo(f"Error running MultiMSOAR: {e}", file=sys.stderr)
        click.echo(f"  Stdout: {e.stdout}", file=sys.stderr)
        click.echo(f"  Stderr: {e.stderr}", file=sys.stderr)
        raise click.Abort()

    # Run summarise_orthogroups_internet_OUT.py
    summarize_script_path = os.path.join(orthonome_dir, "summarise_orthogroups_internet_OUT.py")
    if not os.path.exists(summarize_script_path):
        click.echo(f"Error: Script 'summarise_orthogroups_internet_OUT.py' not found at '{summarize_script_path}'.", file=sys.stderr)
        raise click.Abort()

    genelists_path = os.path.join(work_dir, "genelists.txt")
    spp_list_idx_path = os.path.join(work_dir, "Spp_list.idx")
    output_prefix_path = os.path.join(work_dir, "Orthonome_out") # This is a prefix, script will add extensions.

    summarize_command = [
        sys.executable, summarize_script_path,
        "-l", genelists_path,
        "-m", clusters_file_path, # Allruns.clusters from blast_pairs_dir
        "-i", geneinfo_output_path, # Geneinfo from multi_msoar_inputs_dir
        "-g", orthogroups_output_path, # Orthogroups from multi_msoar_inputs_dir
        "-s", spp_list_idx_path,
        "-o", output_prefix_path
    ]
    click.echo("Running orthogroup summarization...")
    click.echo(f"Command: {' '.join(summarize_command)}")
    try:
        # The original script seems to expect to be run from where the output files should land,
        # or paths should be relative to its CWD. Running from work_dir for simplicity.
        result_summarize = subprocess.run(summarize_command, cwd=work_dir, check=True,
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        click.echo("Orthogroup summarization completed successfully.")
        if result_summarize.stdout: click.echo(f"  Stdout:\n{result_summarize.stdout}")
        if result_summarize.stderr: click.echo(f"  Stderr:\n{result_summarize.stderr}")
    except subprocess.CalledProcessError as e:
        click.echo(f"Error running summarise_orthogroups_internet_OUT.py: {e}", file=sys.stderr)
        click.echo(f"  Stdout: {e.stdout}", file=sys.stderr)
        click.echo(f"  Stderr: {e.stderr}", file=sys.stderr)
        raise click.Abort()
        
    click.echo("run-multimsoar-and-summarize subcommand completed successfully.")

if __name__ == '__main__':
    main_group()
