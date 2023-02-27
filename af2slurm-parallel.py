#!python
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def parse_cmd_args():
    """Sets up argument parser and returns the parsed arguments"""
    parser = ArgumentParser(
        prog="af2slurm-parallel",
        description="Accepts a fasta file and sends an array task to slurm",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", help="Path to fasta file", type=str)
    parser.add_argument("out_dir", help="Directory to write the results to")
    parser.add_argument(
        "--dry-run",
        help="Do not submit any jobs to slurm, just print slurm commands",
        default=False,
        action="store_true",
    )

    ### Control grouping
    parser.add_argument(
        "--max-group-size", help="Put aat most X sequences into one task", type=int, default=30
    )
    parser.add_argument(
        "--max-group-size-AA",
        help="Total number of amino acids in task should be less than X.",
        type=int,
        default=10000,
    )
    parser.add_argument(
        "--max-size-change",
        help="Make a new group if length of amino acids between constructs changes mora than X. Should be synced with recompile_padding",
        type=int,
        default=10,
    )
    #####

    ### Control slurm
    parser.add_argument(
        "--job-name", help="Slurm job name. If empty defaults to fasta name", type=str, default=""
    )
    parser.add_argument(
        "--output",
        help="Slurm output file. If empty defaults to fasta name with .out extension",
        type=str,
        default="",
    )
    parser.add_argument("--partition", help="Slurm partition.", type=str, default="gpu")
    parser.add_argument("--gres", help="GPU to use.", type=str, default="gpu:A40:1")
    #assumes --ntasks=1 is also given
    parser.add_argument("--cpus-per-task", help="How many cpus per task", type=int, default=2)
    

    ### Colabfold batch settings
    
    parser.add_argument("--stop-at-score",
        help="Compute models until plddt (single chain) or ptmscore (complex) > threshold is reached. "
        "This can make colabfold much faster by only running the first model for easy queries.",
        type=float,
        default=100,
    )
    parser.add_argument("--num-recycle",
        help="Number of prediction recycles."
        "Increasing recycles can improve the quality but slows down the prediction.",
        type=int,
        default=None,
    )
    parser.add_argument("--recycle-early-stop-tolerance",
        help="Specify convergence criteria."
        "Run until the distance between recycles is within specified value.",
        type=float,
        default=None,
    )
    parser.add_argument("--num-ensemble",
        help="Number of ensembles."
        "The trunk of the network is run multiple times with different random choices for the MSA cluster centers.",
        type=int,
        default=1,
    )
    parser.add_argument("--num-seeds",
        help="Number of seeds to try. Will iterate from range(random_seed, random_seed+num_seeds)."
        ".",
        type=int,
        default=1,
    )
    parser.add_argument("--random-seed",
        help="Changing the seed for the random number generator can result in different structure predictions.",
        type=int,
        default=0,
    )
    parser.add_argument("--num-models", type=int, default=5, choices=[1, 2, 3, 4, 5])
    parser.add_argument("--recompile-padding",
        type=int,
        default=10,
        help="Whenever the input length changes, the model needs to be recompiled."
        "We pad sequences by specified length, so we can e.g. compute sequence from length 100 to 110 without recompiling."
        "The prediction will become marginally slower for the longer input, "
        "but overall performance increases due to not recompiling. "
        "Set to 0 to disable.",
    )
    parser.add_argument("--model-order", default="1,2,3,4,5", type=str)
    #parser.add_argument("--host-url", default='')
    parser.add_argument("--msa-mode",
        default="mmseqs2_uniref_env",
        choices=[
            "mmseqs2_uniref_env",
            "mmseqs2_uniref",
            "single_sequence",
        ],
        help="Using an a3m file as input overwrites this option",
    )
    parser.add_argument("--model-type",
        help="predict strucutre/complex using the following model."
        'Auto will pick "alphafold2_ptm" for structure predictions and "alphafold2_multimer_v3" for complexes.',
        type=str,
        default="auto",
        choices=[
            "auto",
            "alphafold2",
            "alphafold2_ptm",
            "alphafold2_multimer_v1",
            "alphafold2_multimer_v2",
            "alphafold2_multimer_v3",
        ],
    )
    parser.add_argument("--amber",
        default=False,
        action="store_true",
        help="Use amber for structure refinement."
        "To control number of top ranked structures are relaxed set --num-relax.",
    )
    parser.add_argument("--num-relax",
        help="specify how many of the top ranked structures to relax using amber.",
        type=int,
        default=0,
    )
    parser.add_argument("--templates", default=False, action="store_true", help="Use templates from pdb")
    parser.add_argument("--custom-template-path",
        type=str,
        default=None,
        help="Directory with pdb files to be used as input",
    )
    parser.add_argument("--rank",
        help="rank models by auto, plddt or ptmscore",
        type=str,
        default="auto",
        choices=["auto", "plddt", "ptm", "iptm", "multimer"],
    )
    parser.add_argument("--pair-mode",
        help="rank models by auto, unpaired, paired, unpaired_paired",
        type=str,
        default="unpaired_paired",
        choices=["unpaired", "paired", "unpaired_paired"],
    )
    parser.add_argument("--sort-queries-by",
        help="sort queries by: none, length, random",
        type=str,
        default="length",
        choices=["none", "length", "random"],
    )
    parser.add_argument("--save-single-representations",
        default=False,
        action="store_true",
        help="saves the single representation embeddings of all models",
    )
    parser.add_argument("--save-pair-representations",
        default=False,
        action="store_true",
        help="saves the pair representation embeddings of all models",
    )
    parser.add_argument("--use-dropout",
        default=False,
        action="store_true",
        help="activate dropouts during inference to sample from uncertainty of the models",
    )
    parser.add_argument("--max-seq",
        help="number of sequence clusters to use",
        type=int,
        default=None,
    )
    parser.add_argument("--max-extra-seq",
        help="number of extra sequences to use",
        type=int,
        default=None,
    )
    parser.add_argument("--max-msa",
        help="defines: `max-seq:max-extra-seq` number of sequences to use",
        type=str,
        default=None,
    )
    parser.add_argument("--disable-cluster-profile",
        default=False,
        action="store_true",
        help="EXPERIMENTAL: for multimer models, disable cluster profiles",
    )
    parser.add_argument("--zip",
        default=False,
        action="store_true",
        help="zip all results into one <jobname>.result.zip and delete the original files",
    )
    parser.add_argument("--use-gpu-relax",
        default=False,
        action="store_true",
        help="run amber on GPU instead of CPU",
    )
    parser.add_argument("--save-all",
        default=False,
        action="store_true",
        help="save ALL raw outputs from model to a pickle file",
    )
    parser.add_argument("--save-recycles",
        default=False,
        action="store_true",
        help="save all intermediate predictions at each recycle",
    )
    parser.add_argument("--overwrite-existing-results", default=False, action="store_true")
    parser.add_argument("--disable-unified-memory",
        default=False,
        action="store_true",
        help="if you are getting tensorflow/jax errors it might help to disable this",
    )

    args = parser.parse_args()
    return args

"""
default parameters:
python af2slurm-parallel.py <path/to/fasta/file> <output/directory> \

    --dry-run False \

    ### Control grouping ###
    --max-group-size 30 \
    --max-group-size-AA 10000 \
    --max-size-change 10 \

    ### Slurm controls ###
    --job-name "" \
    --output "" \
    --partition gpu \
    --gres gpu:A40:1 \
    --cpus-per-task 2 \

    ### Colabfold batch settings ###
    --stop-at-score 100 \
    --num-recycle None \
    --recycle-early-stop-tolerance None \
    --num-ensemble 1 \
    --num-seeds 1 \
    --random-seed 0 \
    --num-models 5 \
    --recompile-padding 10 \
    --model-order 1,2,3,4,5 \
    --msa-mode mmseqs2_uniref_env \
    --model-type auto \
    --amber False \
    --num-relax 0 \
    --templates False \
    --custom-template-path None \
    --rank auto \
    --pair-mode unpaired_paired \
    --sort-queries-by length \
    --save-single-representations False \
    --save-pair-representations False \
    --use-dropout False \
    --max-seq None \
    --max-extra-seq None \
    --max-msa None \
    --disable-cluster-profile False \
    --zip False \
    --use-gpu-relax False \
    --save-all False \
    --save-recycles False \
    --overwrite-existing-results False \
    --disable-unified-memory False
"""

def main():
    args = parse_cmd_args()
    print(args)


if __name__ == "__main__":
    main()
