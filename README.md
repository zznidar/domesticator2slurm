# alphafold-slurm-runner
Scripts to run AF2 with slurm integration.


# Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install ConfigArgParse biopython
```

# Usage of alphafold2slurm
TBW


# domesticator-slurm-runner
Scripts to run Domesticator with slurm integration.

# Usage of domesticator2slurm
Copy `.fasta` or `.pdb` file to the specified `in` folder.  
The copied file must start with a **comment** which specifies which **vector**.gb to use and any additional arguments to be passed to Domesticator.

The vector file will be searched for in the specified `vectors` folder. 

## Example:
```fasta
# pET29b.gb
>P1__P2Args
SPEDEIQALEEENAQLEQENAALEEEIAQLEY
```

### Passing additional arguments:
```fasta
# pET29b.gb --nstruct 2 --max_tries 5
>P1__P2Args
SPEDEIQALEEENAQLEQENAALEEEIAQLEY
```