#!python
from configargparse import ArgParser, ArgumentDefaultsHelpFormatter
import subprocess
from time import sleep
from glob import glob
import os
import shutil
from pathlib import Path
import re
import logging
from typing import Tuple


def copy_protein_files(in_path: str, out_folder: str, dry_run: bool = False) -> list:
    """Copies fasta file to the out folder and parses it.
    Returns the path to the parsed copy, output (sub)folder and the first-line comment without the comment char to be used as args for domesticator. This includes the vector filename.extension
    """
    # .fasta || .pdb
    in_path = Path(in_path)
    out_folder = Path(out_folder)

    stem_name = in_path.stem  # stem is file name without extension
    if stem_name.endswith(".fasta"):  # remove the extra .fasta if this is a txt.fasta
        stem_name = stem_name[: -len(".fasta")]

    out_subfolder = out_folder / stem_name
    out_path = out_subfolder / in_path.name

    os.makedirs(out_subfolder, exist_ok=True)

    #copy original input file (with arguments) and add ".original"
    shutil.copy(in_path, str(out_path)+".original")

    # Extract colab_args
    with open(in_path, "r") as source_file:
        # skip empty lines at the start of the file
        lines = source_file.read()
        lines = lines.lstrip(" \n").splitlines()
        if len(lines) == 0:
            logging.warning(f"WARNING: {in_path} is an empty file!")
            return None, None, None
        first_line = lines[0].strip()

        match_colabfold_args_line = re.compile(r"^\s*#\s*")
        if match_colabfold_args_line.match(first_line):  # if the first line matches the pattern
            dom_args = first_line.lstrip(
                "#"
            ).strip()  # Remove only the # symbol from the beginning of the line
            del lines[0]  # Remove the first line from the list
        else:
            dom_args = "" # This should never happen -- vector.gb file is mandatory!
            logging.warning(f"WARNING: {in_path} does not contain # vector.gb! This is not allowed. Please add arguments to the first line of the file.")

    # add fasta header if it is missing. Just use the name of the file. DO NOT DO THIS ON .pdb files
    if lines[0][0] != ">" and in_path.suffix != ".pdb" and in_path.suffix != ".PDB":
        lines.insert(0, ">" + stem_name)

    def filter_stars_spaces(line):
        if len(line)>0 and line[0] == ">":  # if fasta header, don't do any replacements
            return line
        # get rid of stars and spaces in the sequence
        return line.replace("*", "").replace(" ", "")

    if in_path.suffix != ".pdb" and in_path.suffix != ".PDB":  # I don't think .a3m is supported by domesticator?
        lines = [filter_stars_spaces(l) for l in lines]

    # Change out extension to lowercase
    out_path = out_path.with_suffix(out_path.suffix.lower())

    with open(out_path, "w+") as target_file:
        target_file.write("\n".join(lines))
    
    if not dry_run:
        os.remove(in_path)
        
    return out_path, out_subfolder, dom_args

    
def create_slurm_submit_line(protein_path, slurm_options, domesticator_command):
    protein_path = Path(protein_path) # only used to name the job
    # submit line is composed of: sbatch + slurm_args (config file),
    slurm = f"{slurm_options} --parsable --job-name={protein_path.stem} --output={protein_path.with_suffix('.out')} -e {protein_path.with_suffix('.out')} "
    return f"""sbatch  {slurm} --wrap="{domesticator_command}" """

def submit_job(protein_path, args, dom_args, dry_run=False):
    # Vector filename is first argument in dom_args; we precede it with the path to the vectors folder

    dom_command = f"source {args.env_setup_script} && {args.colabfold_path} '{protein_path}' {args.vectors_folder}/{dom_args} --no_idt"

    submit = create_slurm_submit_line(protein_path, args.slurm_args, dom_command)

    if not dry_run:
        slurm_id = subprocess.getoutput(submit)
        logging.info(f"Submitted to slurm with ID {slurm_id}")
    else:
        logging.info(submit)
    

def main():
    parser = ArgParser(
        prog="dom2slurm-watcher",
        description="Watches a folder for gb and fasta files and submits them to slurm",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", help="Path to config file", is_config_file=True, default="dom2slurm.config"
    )
    parser.add_argument(
        "--dry-run",
        # dry-run: program should not actually submit any jobs to a slurm scheduler but should only print out the slurm comand that would be used to submit a job
        help="Do not submit any jobs to slurm, but copy over fasta files and print slurm commands",
        default=False,
        action="store_true",
    )
    parser.add_argument("--in_folder", help="Directory to watch for gb and fasta files", default="./in")
    parser.add_argument("--out_folder", help="Directory to write results to", default="./out")
    parser.add_argument("--log_path_name", help="Filename path to write the log to", default="out.log")
    parser.add_argument("--scan_interval_s", help="Scan folder every X seconds", default=60, type=int)
    parser.add_argument("--vectors_folder", help="Directory with vector.gb files", default="./vectors")
    parser.add_argument(
        "--colabfold_path",
        help="find path to the colabfold",
        default="/home/aljubetic/conda/envs/domesticator/bin/python /home/aljubetic/gits/domesticator3/domesticator3 ",
    )
    parser.add_argument(
        "--slurm_args",
        help="arguments for slurm",
        default="--partition=amd --ntasks=1 --cpus-per-task=1",
    )
    parser.add_argument(
        "--env_setup_script",
        help="python environment that has Domesticator setup",
        default="/home/aljubetic/bin/setup_proxy_settings.sh",
    )
    args = parser.parse_args()

    logging.basicConfig(
        #encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(args.log_path_name), logging.StreamHandler()],
    )

    # Because Domesticator does not support /out folder, change potentially rel /in path to abs path, then change pwd to to the specified out_folder
    args.in_folder = str(Path(args.in_folder).resolve())
    args.out_folder = str(Path(args.out_folder).resolve()) # change out folder as well, in case we refer to it at a later point in time
    args.vectors_folder = str(Path(args.vectors_folder).resolve()) # is without trailing slash
    starting_cwd = os.getcwd()

    logging.info("Running dom2slurm watcher with arguments: " + str(args))

    while True:
        # moves files to the output folder and submit them to slurm
        extensions_prot = [".fasta", ".pdb", ".fasta.txt", ".FASTA", ".PDB"]
        fastas = sorted([f for ext in extensions_prot for f in glob(f"{args.in_folder}/*{ext}")])
        for fasta in fastas:
            logging.info(f'Submitting protein file: "{fasta}"')
            out_protein, out_folder, dom_args = copy_protein_files(fasta, args.out_folder, dry_run=args.dry_run)

            if (out_protein, out_folder, dom_args) == (None, None, None):
                logging.info(f"Skipping {fasta} because it is empty")
                # Rename the empty file to .empty to avoid further processing
                os.rename(fasta, fasta+".empty")    
                continue

            # Wokaround for Domesticator not having --out param 
            os.chdir(out_folder)
            submit_job(out_protein, args, dom_args, dry_run=args.dry_run)
            os.chdir(starting_cwd)

        if args.dry_run:
            # only execute loop once if we are doing a dry run
            break
        sleep(args.scan_interval_s)  # continuous scanning of the input folder for fasta files


if __name__ == "__main__":
    main()
