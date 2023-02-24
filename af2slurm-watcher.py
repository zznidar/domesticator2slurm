#!python
from configargparse import ArgParser, ArgumentDefaultsHelpFormatter
import subprocess
from time import sleep
from glob import glob
import os
import shutil
from pathlib import Path


def copy_over_fasta_file(
    file_path: str,  out_folder: str, dry_run: bool = False
) -> tuple[str, str, str]:  # dry_run is set to 0 which means the comand goes through
    """Makes a folder in the out_folder with the name of the fasta file. Copies fasta file from the in folder the new target folder.
    Returns the path to the fasta file, the output folder and any first line comments, without the comment char (to be used as arguments for colabfold_batch)
    """
    file_path = Path(file_path)
    out_folder = Path(out_folder) 
    if not file_path.suffix == '.fasta':
        raise ValueError("File name does not end with '.fasta'")

    stem_name = file_path.stem # stem is file name without extension 
    
    out_sub_folder = out_folder/stem_name
    out_pathname = out_sub_folder/file_path.name 
    os.makedirs(out_sub_folder, exist_ok=True)
    
    with open(file_path) as source_file:
        #skip empty lines at the start of the file
        lines = source_file.read()
        lines = lines.lstrip(' \n').splitlines()
        first_line = lines[0].strip()
        if first_line[0] == '#': # if the first non empty char of the first line is a comment symbol
            colab_args = first_line[1:]
            del lines[0]
        else:
            colab_args = 0
    
    with open(out_pathname, 'w+') as target_file:
        target_file.write("\n".join(lines))

    return out_pathname, out_sub_folder, colab_args

def create_slurm_submit_line(file_name, slurm_options, colabfold_options):
    # name is just name of fasta file without the .fasta --> for naming
    name = os.path.splitext(file_name)[0]
    # submit line is composed of: sbatch + slurm_args (config file),
    slurm = f"{slurm_options} --job-name={name} --output={name}.out -e {name}.err "
    return f"""sbatch  {slurm} --wrap="{colabfold_options}" """


def move_and_submit_fasta(fasta_path, args, dry_run=False):
    # fasta is a full path to a fasta file in ./in directory
    file_name = os.path.basename(fasta_path)
    target_fasta, out_path_name, colabfold_arguments = copy_over_fasta_file(fasta_path, args.out_folder)

    colabfold_command= f"source {args.env_setup_script} && {args.colabfold_path} {colabfold_arguments} {target_fasta} {out_path_name}"


    submit = create_slurm_submit_line(file_name, args.slurm_args, colabfold_command)
    if not dry_run:
        subprocess.getoutput(submit)
    else:
        print(submit)


def main():
    parser = ArgParser(
        prog="af2slurm-watcher",
        description="Watches a folder for fasta files and submits them to slurm",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", help="Path to config file", is_config_file=True, default="af2slurm.config"
    )
    # config arument specifies path to a configuration file which will be used to set program options
    # is_config_file =True this isn't command line argument(?)
    parser.add_argument(
        "--dry-run",
        # dry-run: program should not actually submit any jobs to a slurm scheduler but should only print out the slurm comand that would be used to submit a job
        help="Do not submit any jobs to slurm, but copy over fasta files and print slurm commands",
        default=True,
        action="store_true",
    )
    parser.add_argument("--in_folder", help="Directory to watch for fasta files", default="./in")
    parser.add_argument("--out_folder", help="Directory to write results to", default="./out")
    parser.add_argument("--scan_interval_s", help="Scan folder every X seconds", default=60, type=int)
    parser.add_argument("--python_path")
    parser.add_argument("--colabfold_path")
    parser.add_argument("--slurm_args")
    parser.add_argument("--env_setup_script")
    args = parser.parse_args()
    # pars_args( method is called on the argument parser 'parser' to actually parse the command-line arguments and populate args variable with the values specified )

    # print(args)

    while True:
        # inf loop in which it scans the input directory for new fasta files and calls one of the functions from above - move and submit fasta
        # moves files to the output folder and submit them to slurm
        fastas = sorted(glob(f"{args.in_folder}/*.fasta"))
        for fasta in fastas:
            print(f"Submitting file: {fasta}")
            move_and_submit_fasta(fasta, args, dry_run=args.dry_run)
        if args.dry_run:
            # only execute loop once if we are doing a dry run
            break
        sleep(args.scan_interval_s)  # continuous scanning of the input folder for fasta files


if __name__ == "__main__":
    main()
