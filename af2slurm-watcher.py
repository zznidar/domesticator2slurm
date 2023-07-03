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


def move_over_fasta_file(
    file_path: str, out_folder: str, dry_run: bool = False
) -> Tuple[str, str, str]:  # dry_run is set to 0 which means the comand goes through
    """Makes a file to the out folder.
    Returns the path to the fasta/MSA file, the output folder and any first line comments, without the comment char (to be used as arguments for colabfold_batch)
    """
    file_path = Path(file_path)
    out_folder = Path(out_folder)

    stem_name = file_path.stem  # stem is file name without extension
    if stem_name.endswith(".fasta"):  # remove the extra .fasta if this is a txt.fasta
        stem_name = stem_name[: -len(".fasta")]

    out_sub_folder = out_folder / stem_name
    out_pathname = out_sub_folder / file_path.name
    os.makedirs(out_sub_folder, exist_ok=True)

    #copy original input file (with arguments) and add ".original"
    shutil.copy(file_path, str(out_pathname)+".original")

    # Extract colab_args
    with open(file_path) as source_file:
        # skip empty lines at the start of the file
        lines = source_file.read()
        lines = lines.lstrip(" \n").splitlines()
        first_line = lines[0].strip()

        match_colabfold_args_line = re.compile(r"^\s*#\s*-\s*")
        if match_colabfold_args_line.match(first_line):  # if the first line matches the pattern
            colab_args = first_line.lstrip(
                "#"
            ).strip()  # Remove only the # symbol from the beginning of the line
            del lines[0]  # Remove the first line from the list
        else:
            colab_args = ""

    # add fasta header if it is missing. Just use the name of the file. DO NOT DO THIS ON .a3m files
    if lines[0][0] != ">" and file_path.suffix != ".a3m":
        lines.insert(0, ">" + stem_name)

    def filter_stars_spaces(line):
        if len(line)>0 and line[0] == ">":  # if fasta header, don't do any replacements
            return line
        # get rid of stars and spaces in the sequence
        return line.replace("*", "").replace(" ", "")

    if file_path.suffix != ".a3m":  # no filtering should be needed for .a3m files
        lines = [filter_stars_spaces(l) for l in lines]

    with open(out_pathname, "w+") as target_file:
        target_file.write("\n".join(lines))

    if not dry_run:
        os.remove(file_path)

    return out_pathname, out_sub_folder, colab_args


def create_slurm_submit_line(file_name, slurm_options, colabfold_options):
    file_name = Path(file_name)
    # submit line is composed of: sbatch + slurm_args (config file),
    slurm = f"{slurm_options} --parsable --job-name={file_name.stem} --output={file_name.with_suffix('.out')} -e {file_name.with_suffix('.out')} "
    return f"""sbatch  {slurm} --wrap="{colabfold_options}" """


def move_and_submit_fasta(fasta_path, args, dry_run=False):
    # fast_path is a full path to a fasta file in ./in directory
    target_fasta, out_path_name, colabfold_arguments = move_over_fasta_file(
        fasta_path, out_folder=args.out_folder, dry_run=args.dry_run
    )

    colabfold_command = f"source {args.env_setup_script} && {args.colabfold_path} {colabfold_arguments} {target_fasta} {out_path_name}"

    submit = create_slurm_submit_line(target_fasta, args.slurm_args, colabfold_command)

    if not dry_run:
        slurm_id = subprocess.getoutput(submit)
        logging.info(f"Submitted to slurm with ID {slurm_id}")
    else:
        logging.info(submit)


def main():
    parser = ArgParser(
        prog="af2slurm-watcher",
        description="Watches a folder for fasta files and submits them to slurm",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", help="Path to config file", is_config_file=True, default="af2slurm.config"
    )
    parser.add_argument(
        "--dry-run",
        # dry-run: program should not actually submit any jobs to a slurm scheduler but should only print out the slurm comand that would be used to submit a job
        help="Do not submit any jobs to slurm, but copy over fasta files and print slurm commands",
        default=False,
        action="store_true",
    )
    parser.add_argument("--in_folder", help="Directory to watch for fasta files", default="./in")
    parser.add_argument("--out_folder", help="Directory to write results to", default="./out")
    parser.add_argument("--log_path_name", help="Directory to write results to", default="out.log")
    parser.add_argument("--scan_interval_s", help="Scan folder every X seconds", default=60, type=int)
    parser.add_argument("--nochange_interval_s", help="When scanning, wait X seconds to confirm all files have been copied", default=5, type=int)
    parser.add_argument(
        "--colabfold_path",
        help="find path to the colabfold",
        default="/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch ",
    )
    parser.add_argument(
        "--slurm_args",
        help="arguments for slurm",
        default="--partition=gpu --gres=gpu:A40:1 --ntasks=1 --cpus-per-task=2",
    )
    parser.add_argument(
        "--env_setup_script",
        help="python environment that has AF2 setup",
        default="/home/aljubetic/bin/set_up_AF2.3.sh",
    )
    args = parser.parse_args()

    logging.basicConfig(
        #encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(args.log_path_name), logging.StreamHandler()],
    )

    logging.info("Running af2slurm watcher with arguments: " + str(args))

    while True:
        # moves files to the output folder and submit them to slurm
        extensions = [".fasta", ".a3m", ".fasta.txt"]
        fastas_old = []
        while True:
            fastas = sorted([f for ext in extensions for f in glob(f"{args.in_folder}/*{ext}")])
            if fastas == fastas_old:
                logging.info(f"Debounced. Found {len(fastas)} fasta files.") if len(fastas) > 0 else None
                break
            fastas_old = fastas # We can do this because we never append -- we always set the array to a new value (pointer)
            sleep(args.nochange_interval_s)
        for fasta in fastas:
            logging.info(f"Submitting file: {fasta}")
            move_and_submit_fasta(fasta, args, dry_run=args.dry_run)
        if args.dry_run:
            # only execute loop once if we are doing a dry run
            break
        sleep(args.scan_interval_s)  # continuous scanning of the input folder for fasta files


if __name__ == "__main__":
    main()
