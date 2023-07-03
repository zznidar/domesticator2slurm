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

print("Domesticator 3 slurm watcher")

# what we do is copy files, schedule run for each gb, then delete files on our end. 
# TODO: On server end, upon finishing Domesticator run, delete the gb in question, then check the folder for other gbs. Should there be more, do nothing. If there are no gbs left, delete fastas as well.
# Additional arguments should be written in the .gb file since we run Domesticator for each gb (and pass all fastas to it)

def copy_protein_files(in_paths: list, out_folder: str, dry_run: bool = False) -> list:
    """Copies fasta files to the out folder.
    Returns the paths to the copies of fasta files.
    """
    # .fasta || .pdb
    os.makedirs(out_folder, exist_ok=True)
    out_paths = []
    for in_path in in_paths:
        in_path = Path(in_path)
        out_path = Path(out_folder) / in_path.name
        out_paths.append(str(out_path))
        shutil.copy(in_path, out_path)
        if not dry_run:
            os.remove(in_path)
    # TODO: add fasta header if missing!
    print(f"copying prots {in_paths} to {out_paths}")
    return out_paths
    
def copy_vector_file(in_path: str, out_folder: str, dry_run: bool = False) -> Tuple [str, str, str]:
    """Copies gb file to the out folder.
    Returns the path to the gb file, output folder and any first line comments, without the comment char (to be used as arguments for domesticator_batch)
    """
    # .gb
    in_path = Path(in_path)
    out_folder = Path(out_folder)

    stem_name = in_path.stem  # stem is file name without extension
    if stem_name.endswith(".gb"):  # remove the extra .gb if this is a txt.gb
        stem_name = stem_name[: -len(".gb")]

    out_subfolder = out_folder / stem_name
    out_path = out_subfolder / in_path.name

    os.makedirs(out_subfolder, exist_ok=True)

    # Copy original input gb (with args) as *.original
    shutil.copy(in_path, str(out_path) + ".original")
    print(f"copying vec {in_path} to {out_path}")

    # Extract first line comments
    with open(in_path, "r") as source_file:
        # skip empty lines at the beginning of the file
        lines = source_file.read()
        lines = lines.lstrip(" \n").splitlines()
        first_line = lines[0].strip()

        match_domesticator_args_line = re.compile(r"^\s*#\s*-\s*")
        if match_domesticator_args_line.match(first_line):  # if the first line matches the pattern
            dom_args = first_line.lstrip(
                "#"
            ).strip()  # Remove only the # symbol from the beginning of the line
            del lines[0]  # Remove the first line from the list
        else:
            dom_args = ""
        
        with open(out_path, "w+") as target_file:
            target_file.write("\n".join(lines))
        
        if not dry_run:
            os.remove(in_path)
        
    return out_path, out_subfolder, dom_args


def create_slurm_submit_line(vector_path, slurm_options, domesticator_command):
    vector_path = Path(vector_path) # only used to name the job
    # submit line is composed of: sbatch + slurm_args (config file),
    slurm = f"{slurm_options} --parsable --job-name={vector_path.stem} --output={vector_path.with_suffix('.out')} -e {vector_path.with_suffix('.out')} "
    return f"""sbatch  {slurm} --wrap="{domesticator_command}" """

def submit_job(vector_path, protein_paths, args, dom_args, dry_run=False):
    # vector_path is a full path to a gb file in ./in directory

    dom_command = f"source {args.env_setup_script} && {args.colabfold_path} {' '.join(protein_paths)} {vector_path} {dom_args} --no_idt"

    submit = create_slurm_submit_line(vector_path, args.slurm_args, dom_command)

    if not dry_run:
        slurm_id = subprocess.getoutput(submit)
        logging.info(f"Submitted to slurm with ID {slurm_id}")
    else:
        logging.info(submit)
    

## Delete this function, it serves no purpose
def now_delete_fasta(fasta, dry_run=False):
    if dry_run:
        logging.info(f"Would delete {fasta}, but skipping since dry_run is set")
    else:
        file_path = Path(fasta)
        os.remove(file_path)


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
    parser.add_argument("--nochange_interval_s", help="When scanning, wait X seconds to confirm all files have been copied", default=5, type=int)
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
        help="python environment that has AF2 setup",
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
    starting_cwd = os.getcwd()

    # Also, for each vector, we should create a new subfolder -- just because the output is always named order.dna.fasta, which means they get overwritten if multiple vectors are submitted.


    logging.info("Running dom2slurm watcher with arguments: " + str(args))

    while True:
        # moves files to the output folder and submit them to slurm
        extensions_prot = [".fasta", ".pdb", ".fasta.txt"]
        extensions_vec = [".gb"]

        # Debouncer: wait until no new files are added to the folder for X seconds
        files_old = []
        while True:
            fastas = sorted([f for ext in extensions_prot for f in glob(f"{args.in_folder}/*{ext}")])
            gbs = sorted([f for ext in extensions_vec for f in glob(f"{args.in_folder}/*{ext}")])
            if files_old == fastas + gbs:
                logging.info(f"Debounced. Found {len(fastas)} fasta files and {len(gbs)} gb files.")
                break
            files_old = fastas + gbs
            logging.info(f"They are not the same. {files_old} != {fastas + gbs}. Waiting {args.nochange_interval_s} seconds.")
            sleep(args.nochange_interval_s)

        print(f"Found the following files: {fastas} and {gbs}")

        # Copy all proteins. 
        out_proteins = copy_protein_files(fastas, args.out_folder, dry_run=args.dry_run)

        # For each vector, submit all fastas (multiple _protein_ args can be specified, and one _vector_)
        for gb in gbs:
            logging.info(f"Submitting vector file: {gb}")
            out_vector, out_folder, dom_args = copy_vector_file(gb, args.out_folder, dry_run=args.dry_run)

            os.chdir(out_folder)
            submit_job(out_vector, out_proteins, args, dom_args, dry_run=args.dry_run)
            os.chdir(starting_cwd)

        if args.dry_run:
            # only execute loop once if we are doing a dry run
            break
        sleep(args.scan_interval_s)  # continuous scanning of the input folder for fasta files

        # Also add an argument --nochange_interval_s
        # After sleeping, scan for files. Wait --nochange seconds, then rescan. Should files match, copy them. Else: wait another --nochange seconds.


if __name__ == "__main__":
    main()
