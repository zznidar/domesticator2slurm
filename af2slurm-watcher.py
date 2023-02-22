#!python
from configparser import ConfigParser
from argparse import ArgumentParser
from time import sleep
from glob import glob

def copy_over_fasta_file(file_name: str, in_folder: str, out_folder: str, dry_run: bool = True) -> tuple[str, str]:
    """Makes a folder in the out_foder with the name of the fasta fie. Copies fasta file from the in folder the new target folder. Returns the path to the fasta file and any first line comments (to be used as arguments for colabfold_batch)"""
    pass

def create_slurm_submit_line(fast_file, slurm_options, colabfold_options):
    pass

def move_and_submit_fasta(fasta, cf, dry_run=False):
    ""
    #copy_over_fasta_file(fasta)
    pass

def main():
    parser = ArgumentParser(
        prog="af2slurm-watcher", description="Watches a folder for fasta files and submits them to slurm"
    )
    parser.add_argument("--config", help="Path to config file", type=str, default="af2slurm.config")
    parser.add_argument(
        "--dry-run",
        help="Do not submit any jobs to slurm, but copy over fasta files and print slurm commands",
        default=False,
        action="store_true",
    )
    args = parser.parse_args()

    # print(args.config)
    cf = ConfigParser()
    cf.read(args.config)
    # print(cf.get("APP", "in_folder"))

    scan_interval_s=cf.getint('APP', 'scan_interval_s')
    in_folder=cf.get('APP', 'in_folder')

    while True:
        fastas = sorted(glob(f'{in_folder}/*.fasta'))
        for fasta in fastas:
            print(f"Submitting file: {fasta}" )
            move_and_submit_fasta(fasta, cf, dry_run=args.dry_run)
        sleep(scan_interval_s)

    # TODO event loop here
    # until true do scan_foders, sleep

if __name__ == "__main__":
    main()
