## PEP and Looper Config to download reference sequence assets to brickyard

### Uses Looper>2.0.0


#### For each PEP,

- check input_pep, i.e. `pep_config`
- check for the correct output pep so as to not overwrite results, i.e.
```
pipestat:
  pephub_path: donaldcampbelljr/human_seqcol_digests:default
```

#### Current looper_config files:

- .looper.yaml -> for running main set of initially curated files (human)
- .looper_ncbi_38.yaml -> only hg38 samples acquired from NCBI
- .looper_mm.yaml -> initial mouse set
- .looper_local.yaml -> for files that could not be downloaded with the original pipeline but were instead manually downloaded (e.g. via web browser on HPC) and then processed (human)
- .looper_ensembl -> only Ensembl hg38 samples

#### Example commands to run on HPC:

##### load miniforge and conda env with require packages
```
module load miniforge
conda activate seqcolwork
```

##### set looper results directory
```
export RESULTS=/home/zzz3fh/LOOPERRESULTS/seqcol06may2025/local_fasta_to_digests
```

##### login to PEPhub so that results can be reported to a PEP
```
phc login
```

##### Run a looper config
```
looper run -c .looper_local.yaml --package slurm --compute PARTITION=standard time='01-00:00:00' cores='16' mem='16000' -d
```