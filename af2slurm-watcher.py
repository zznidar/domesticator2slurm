#!python
from configparser import ConfigParser
from argparse import ArgumentParser


def copy_over_fasta_file(file_name: str, in_folder: str, out_folder: str) -> tuple[str, str]:
    """Makes a folder in the out_foder with the name of the fasta fie. Copies fasta file from the in folder the new target folder. Returns the path to the fasta file and any first line comments (to be used as arguments for colabfold_batch)"""
    pass

def create_slurm_submit_line(fast_file, slurm_options, colabfold_options):
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

    # TODO event loop here
    # until true do scan_foders, sleep

if __name__ == "__main__":
    main()
