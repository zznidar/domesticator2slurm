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
from datetime import datetime


def move_over_yaml_file(
    file_path: str, out_folder: str, dry_run: bool = False
) -> Tuple[str, str]: #**now only two return values # dry_run is set to 0 which means the comand goes through
    """Copy the yaml file to the out folder.
    Returns the path to the yaml file and the output folder
    """
    file_path = Path(file_path)
    out_folder = Path(out_folder)

    stem_name = file_path.stem  # stem is file name without extension

    out_sub_folder = out_folder / f"{stem_name}_{datetime.now().strftime('%H_%M_%S')}"
    #out_sub_folder = out_folder / stem_name
    out_pathname = out_sub_folder / file_path.name
    # Check if the output subfolder already exists; if so, add "X" to the out_pathname and out_sub_folder
    
    """
    if out_sub_folder.exists():
        out_sub_folder = Path(str(out_sub_folder) + "_X")
        out_pathname = out_sub_folder / file_path.name
"""
    os.makedirs(out_sub_folder, exist_ok=True)

    #copy original input file (with arguments) and add ".original"
    #shutil.copy(file_path, str(out_pathname)+".original")
    
    shutil.copy(file_path, str(out_pathname))

    if not dry_run:
        os.remove(file_path)

    return out_pathname, out_sub_folder

# CHANGE TO PROSCULPT COMMAND -> "sbatch launch.sh run_details.yaml" 
def create_slurm_submit_line(yaml_file_name, launch_script_path,out_path_name, env_setup_script):
    # Submit the file we actually created: the .original copy, with absolute path
    yaml_file_name = Path(yaml_file_name).resolve()
    return f"bash -c 'source {env_setup_script} $$ /home/folivieri/miniforge3/etc/profile.d/conda.sh && conda activate /home/folivieri/.conda/envs/prosculpt && python {launch_script_path} {yaml_file_name} ++output_dir={out_path_name}'"
#TODO: ne rabis env providat


def move_and_submit_yaml(yaml_path, args, dry_run=False):
    # yaml_path is a full path to the yaml file in ./in directory
    target_yaml, out_path_name = move_over_yaml_file(
        yaml_path, out_folder=args.out_folder, dry_run=args.dry_run
    )
    
    submit = create_slurm_submit_line(target_yaml, args.launch_script_path, out_path_name, args.env_setup_script) 

    if not dry_run:
        os.chdir(args.out_folder)  # change to the output directory
        slurm_id = subprocess.getoutput(submit)
        logging.info(f"Submitted to slurm with ID {slurm_id}")
        os.chdir("../../")  # go back to the original directory
    else:
        logging.info(submit)


def main():
    parser = ArgParser(
        prog="ps2slurm-watcher",
        description="Watches a folder for yaml files and submits them to slurm",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", help="Path to config file", is_config_file=True, default="ps2slurm.config"
    )
    parser.add_argument(
        "--dry-run",
        # dry-run: program should not actually submit any jobs to a slurm scheduler but should only print out the slurm comand that would be used to submit a job
        help="Do not submit any jobs to slurm, but copy over yaml files and print slurm commands",
        default=False,
        action="store_true",
    )
    parser.add_argument("--in_folder", help="Directory to watch for yaml files", default="./in")
    parser.add_argument("--out_folder", help="Directory to write results to", default="./out")
    parser.add_argument("--log_path_name", help="Directory to write results to", default="out.log")
    parser.add_argument("--scan_interval_s", help="Scan folder every X seconds", default=60, type=int)
    parser.add_argument(
        "--launch_script_path",
        help="find path to the ps launch script",
        default="/home/folivieri/prosculpt/slurm_runner.py",
    )
    
    parser.add_argument(
        "--env_setup_script",
        help="script to set up prosculpt env and proxy",
        default="/home/d12/ps2slurm/ps2slurm_v1/set_up_prosculpt_env.sh",
    )
    
    args = parser.parse_args()

    logging.basicConfig(
        #encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
        handlers=[logging.FileHandler(args.log_path_name), logging.StreamHandler()],
    )

    logging.info("Running ps2slurm watcher with arguments: " + str(args))
    
    while True:
        # moves files to the output folder and submit them to slurm
        extensions = [".yaml", ".yml", ".YAML"] #** Monitor for yaml files
        yamls = sorted([f for ext in extensions for f in glob(f"{args.in_folder}/*{ext}")])
        logging.info(f"Found YAML files: {yamls}")
        for yaml in yamls:
            logging.info(f"Submitting file: {yaml}")
            move_and_submit_yaml(yaml, args, dry_run=args.dry_run)
        if args.dry_run:
            # only execute loop once if we are doing a dry run
            break
        sleep(args.scan_interval_s)  # continuous scanning of the input folder for files


if __name__ == "__main__":
    main()
