#!python
from configparser import ConfigParser
from argparse import ArgumentParser


if __name__ == "__main__":
    parser = ArgumentParser(
        prog="af2slurm-watcher", description="Watches a folder for fasta files and submits them to slurm"
    )
    parser.add_argument("--config", help="Path to config file", type=str, default="af2slurm.config")
    args = parser.parse_args()

    #print(args.config)
    cf = ConfigParser()
    cf.read(args.config)
    #print(cf.get("APP", "in_folder"))
