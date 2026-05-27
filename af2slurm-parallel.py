#!python
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
import socket
from glob import glob
from pathlib import Path
from Bio import SeqIO
from math import ceil


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
        "--proteinmpnn",
        help="Fasta input is ProteinMPNN output file. Ranks sequences by score",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--filter-proteinmpnn", help="Filter best X fasta sequences sorted by 'Score'", type=int, default=1000
    )
    parser.add_argument(
        "--mpnn-include-original", help="Include original sequence in the fasta", type=bool, default=False
    )
    parser.add_argument(
        "--target",
        help="target sequence to predict along with your designed binder",
        type=str,
        default=None,
    )
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
        "--job-name", help="Slurm job name. If empty defaults to fasta name", type=str, default="af2slurm-parallel"
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
    parser.add_argument(
        "--num-recycle",
        help="Number of prediction recycles. "
        "Increasing recycles can improve the prediction quality but slows down the prediction.",
        type=int,
        default=None,
    )
    parser.add_argument("--recycle-early-stop-tolerance",
        help="Specify convergence criteria."
        "Run until the distance between recycles is within specified value.",
        type=float,
        default=100,
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
    parser.add_argument(
        "--msa-mode",
        default="mmseqs2_uniref_env",
        choices=[
            "mmseqs2_uniref_env",
            "mmseqs2_uniref",
            "single_sequence",
        ],
        help="Databases to use to create the MSA: UniRef30+Environmental (default), UniRef30 only or None. "
        "Using an A3M file as input overwrites this option.",
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

# Define a function to write a list of sequences to a FASTA file
def write_to_fasta(name, list_of_seq):
    with open(name, 'w') as out_file:
        for seq in list_of_seq:
            out_file.write(f'>{seq.id}\n')
            out_file.write(f'{seq.seq}\n')

# Define a function to parse the description line of a FASTA record
def parse_fasta_description(desc):
    parts = [p.strip() for p in desc.split(',')]
    data = {k.strip(): v.strip() for k, v in [p.split('=') for p in parts]}
    return data


def main():
    args = parse_cmd_args()
    
    # Main arguments
    out_dir = Path(args.out_dir)
    fasta_file = args.fasta
    job_name = args.job_name if args.job_name else fasta_file
    proteinmppn = args.proteinmpnn
    target_sequence = args.target
    filter_proteinmpnn = args.filter_proteinmpnn

    # Set maximum group size, total amino acid count, and size change for clustering sequences
    MAX_GROUP_SIZE = args.max_group_size #--max-group-size 30 {args.save_recycles}
    MAX_GROUP_TOTAL_AA = args.max_group_size_AA #--max-group-size-AA
    MAX_SIZE_CHANGE = args.max_size_change #--max-size-change
    GROUPS = []

    # Create the output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    # If ProteinMPNN is input, sort by score and take best X
    if proteinmppn:
         # Read in a list of sequences from a FASTA file and sort them by score
        seqs = list(SeqIO.parse(open(fasta_file), 'fasta'))
        
        # Sort without first sequence
        first = seqs[0]
        seqs = seqs[1:]
        len_before_deduplication = len(seqs)
        #remove duplicates
        seqs_dict = {seq.seq: seq for seq in seqs}
        seqs = list(seqs_dict.values())
        
        len_after_deduplication = len(seqs)
        print(f"Removed {len_before_deduplication-len_after_deduplication} duplicate sequences")
        seqs = sorted(seqs, key=lambda seq: float(parse_fasta_description(seq.description)['score']))
        seqs = [first] + seqs
    
        # Take the top X sequences and the native sequence
        seqs = seqs[:filter_proteinmpnn+1]

        # Write the selected sequences to a new FASTA file with modified IDs
        with open(f'{out_dir}/{job_name}.fasta', 'w') as seq_list:
            if args.mpnn_include_original:
                seq_list.write(f">Original_sequence\n{str(seqs[0].seq).replace(' ','').replace('-','').replace('/',':')}{f':{target_sequence}' if target_sequence else ''}\n")
            for seq in seqs[1:]:
                info = parse_fasta_description(seq.description)
                new_id = f"{info['sample']}|{info['score']} {info['T']} {info['score']}"
                seq_list.write(f">{new_id}\n{str(seq.seq).replace(' ','').replace('-','').replace('/',':')}{f':{target_sequence}' if target_sequence else ''}\n")
        
        # Read in a list of sequences from the new FASTA file
        seq_list = list(SeqIO.parse(open(seq_list.name), 'fasta'))
    else:
        seq_list = list(SeqIO.parse(open(fasta_file), 'fasta'))

    # Initialize a new group with the first sequence
    seq = seq_list[0]
    group_index = 0
    GROUPS.append([seq])
    group_size = 1
    group_total_AA = len(seq)
    last_added_length = len(seq)

    # Iterate over the remaining sequences and cluster them into groups
    for seq in seq_list[1:]:
        # If there is a change in criteria, create a new group
        size_change = len(seq) / last_added_length
        if group_size > MAX_GROUP_SIZE or group_total_AA > MAX_GROUP_TOTAL_AA or size_change > MAX_SIZE_CHANGE:
            group_index += 1
            GROUPS.append([seq])
            group_size = 1
            group_total_AA = len(seq)
            last_added_length = len(seq)
        else:
            GROUPS[group_index].append(seq)
            group_size += 1
            group_total_AA += len(seq)
            last_added_length = len(seq)

    #Write each group to a separate file and generate commands to run the ColabFold program on each file
    for n, g in enumerate(GROUPS):
        write_to_fasta(out_dir / f"g{n:04d}.fasta", g)

    fastas = sorted(glob(f'{out_dir}/g*.fasta'))

    #Create a file to store the commands
    #print(f'{out_dir}/run.tasks')
    #TODO confing should be read from config file and not hardcoded
    with open(f'{out_dir}/run.tasks', 'w') as f:
        for fasta in fastas:
            fasta_name = Path(fasta).stem
            f.write(f'. /home/aljubetic/bin/set_up_AF2.3.sh && mkdir -p {out_dir / fasta_name} && ' \
                f'/home/aljubetic/AF2/CF2.3/colabfold-conda/bin/colabfold_batch ' \
                #f'--stop-at-score {args.stop_at_score} ' \
                #f'--num-recycle {args.num_recycle} ' \
                #f'--recycle-early-stop-tolerance {args.recycle_early_stop_tolerance} ' \
               # f'--num-ensemble {args.num_ensemble} ' \
               # f'--num-seeds {args.num_seeds} ' \
               # f'--random-seed {args.random_seed} ' \
               # f'--num-models {args.num_models} ' \
               # f'--recompile-padding {args.recompile_padding} ' \
               # f'--model-order {args.model_order} ' \
                f'--msa-mode {args.msa_mode} ' \
                f'--model-type {args.model_type} ' \
               # f'--amber {args.amber} ' \
               # f'--num-relax {args.num_relax} ' \
               # f'--templates {args.templates} ' \
               # f'--custom-template-path {args.custom_template_path} ' \
               # f'--rank {args.rank} ' \
               # f'--pair-mode {args.pair_mode} ' \
               # f'--sort-queries-by {args.sort_queries_by} ' \
               # f'--save-single-representations {args.save_single_representations} ' \
               # f'--save-pair-representations {args.save_pair_representations} ' \
               # f'--use-dropout {args.use_dropout} ' \
               # f'--max-seq {args.max_seq} ' \
               # f'--max-extra-seq {args.max_extra_seq} ' \
               # f'--max-msa {args.max_msa} ' \
               # f'--disable-cluster-profile {args.disable_cluster_profile} ' \
               # f'--zip {args.zip} ' \
               # f'--use-gpu-relax {args.use_gpu_relax} ' \
               # f'--save-all {args.save_all} ' \
               # f'--save-recycles {args.save_recycles} ' \
               # f'--overwrite-existing-results {args.overwrite_existing_results} ' \
               # f'--disable-unified-memory {args.disable_unified_memory}' \
               f'{fasta} {out_dir / fasta_name} ' \
                '\n')

    #Read the commands from the file
    with open(f'{out_dir}/run.tasks') as cmds:
        lines = cmds.read()
        lines = lines.split('\n')
    # delete empty lines
    lines = [line for line in lines if line.strip()]
    # Prepare params for the jobs
    job_name = args.job_name if args.job_name else fasta_file
    output_file = args.output if args.output else f"{job_name}.out"
    slurm_params =  f'--partition={args.partition} --gres={args.gres} --ntasks=1 ' \
                    f'--cpus-per-task={args.cpus_per_task} --job-name={job_name} ' \
                    f'--output={out_dir}/{output_file} -e {out_dir}/{job_name}.err '

    GROUP_SIZE=1
    task_list = f'{out_dir}/run.tasks'
    num_tasks = ceil(len(lines)/GROUP_SIZE)
    print(num_tasks)
    
    #Submit
    dry_run = args.dry_run
    

    cmd_string = f"export GROUP_SIZE=1; sbatch {slurm_params} -a 1-{num_tasks} /home/aljubetic/scripts/wrapper_slurm_array_job_group.sh {task_list}"
    if dry_run:
        print(cmd_string)
    else:
        os.system(cmd_string)


if __name__ == "__main__":
    main()