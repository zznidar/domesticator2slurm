#!python
from configargparse import ArgParser, ArgumentDefaultsHelpFormatter
import subprocess
from time import sleep
from glob import glob
import os
import shutil


def copy_over_fasta_file(
    file_name: str, in_folder: str, out_folder: str, dry_run: bool = False
) -> tuple[str, str]:  # dry_run is seto to 0 which means the comand goes through
    """Makes a folder in the out_foder with the name of the fasta file. Copies fasta file from the in folder the new target folder. Returns the path to the fasta file and any first line comments (to be used as arguments for colabfold_batch)"""
    # Check if the file name ends with '.fasta'
    if not file_name.endswith(".fasta"):
        raise ValueError("File name does not end with '.fasta'")

    # Create input path by joining in_folder (from arguments in config file) and filename with ending .fasta
    in_path = os.path.join(in_folder, file_name)
    # create outpath y joining out_folder (from arguments file) and filename but without ending .fasta
    out_path = os.path.join(out_folder, os.path.splitext(file_name)[0])
    # if out_path doesn't exist yet make a folder in the ./out directory (out_path) and copy over the fasta file
    if not os.path.exists(out_path):
        os.makedirs(out_path, exist_ok=True)
        print("creating folder")
        if not dry_run:
            shutil.copy(in_path, out_path)
        # after copying the fasta into the ./out/name/ direcotry open it and read the first line
        with open(in_path, "r") as f:
            first_line = f.readline().strip()
            # if the line has a first character different than '>' save it to the args variable, if the first lien starts with '>' than args variable is empty
            if first_line[0] != ">":
                args = first_line[1:]
            else:
                args = "no arguments"

        # Return two variables args (contains first line arguments from the fasta file) and out_path ./out/file_name
        return (args, out_path)
    # if the path already exists returns two empty variables     !!!!optimize this part!!!!!! - idalno da preveri že na začetku kateri .fasta fili so novo ali pa da jih sproti, ko jih kopira tudi briše iz ./in folderja
    return ("", "")


def create_slurm_submit_line(file_name, slurm_options, colabfold_options):
    # name is just name of fasta file without the .fasta --> for naming
    name = os.path.splitext(file_name)[0]
    # submit line is composed of: sbatch + slurm_args (config file),
    slurm = f"{slurm_options} --job-name={name} --output={name}.out -e {name}.err "
    return f"""sbatch  {slurm} --wrap="{colabfold_options}" """

    pass


def move_and_submit_fasta(fasta, args, dry_run=False):
    # fasta is a full path to a fasta file in ./in directory
    file_name = os.path.basename(fasta)
    first_line, out_path = copy_over_fasta_file(file_name, args.in_folder, args.out_folder)
    # first line is currently not included in the AF2 comand, need to check why --save_recycles gets marked as invalid flag for AF2!!!!! (ko to ugotovim - dodaj {first_line} za pythn_path in verjetno izbriši --msa-mode)
    # create a string for colabfold options NEED SOME MORE WORK
    colabfold_options = f"source {args.env_setup_script} && {args.python_path} --msa-mode single_sequence {out_path}/{file_name} {out_path}"
    #
    if out_path != "":
        submit = create_slurm_submit_line(file_name, args.slurm_args, colabfold_options)
        subprocess.getoutput(submit)

    ""
    # copy_over_fasta_file(fasta)
    pass


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

    for i in range(3):
        # while True:
        # inf loop in which it scans the input directory for new fasta files and calls one of the functions from above - move and submit fasta
        # moves files to the output folder and submit them to slurm
        fastas = sorted(glob(f"{args.in_folder}/*.fasta"))
        for fasta in fastas:
            print(f"Submitting file: {fasta}")
            move_and_submit_fasta(fasta, args, dry_run=args.dry_run)
        sleep(args.scan_interval_s)  # continous scanning of the input foldrer for fasta files

    # TODO event loop here
    # until true do scan_foders, sleep


if __name__ == "__main__":
    main()
