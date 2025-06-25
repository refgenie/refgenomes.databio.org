# Sequence Collection Comparison Pipeline

### Sample Level Pipeline

For each sample in your PEP:
- create sequence_collection dictionary and dump it to a json
- create chrom.sizes files from the sequence_collection dict

### Project Level Pipeline

Calculate stats for all samples in the PEP including, jaccard similarities for names, lengths, name_len_pairs, sequences.

### Plotting statistics figures
After running project pipeline and calculating statistics, you can plot figures. 

You need to provide the PEPhub path for the PEP with stats results, species (for figure titles), PEPhub path that links sample_name with authority, path to local directory with jsons of sequence collection information.

You must also have a folder of json representations of each of the reference genomes you wish to analyze and plot (these jsons are created via the above sample-level pipeline).

Example command
```
python3 stats_graphs.py /home/drc/Downloads/refgenomes_pics_test/22may2025/mus_musculus/ donaldcampbelljr/mouse_seq_col_results:default mus_musculus /home/drc/Downloads/mouse_jsons_from_rivanna/json/ donaldcampbelljr/mouse_seqcol_digests:default 
```

There is also an R script, `test_plotting.R`. You must add the path to the output folder for the above python script (it produces csvs for plotting)

#### plotting clustered sequence frequency
- run script `cluster_seqcounts.py` in the clustering_scripts folder and then point `cluster_sequence_plots.R` to its output and run to produce a figure of sequences clustered by frequency.


### ENV Variables to Set before execution

RESULTS -> where looper will put its results

BRICKYARD_JSON_LOCATION -> where to place created JSON files

BRICKYARD_CHROM_SIZES_LOCATION -> where to place created chrom.sizes files


## Various PEPs for digests and results

### Homo sapiens

- all sampled Genomes for Homo sapiens (~52)
	- donaldcampbelljr/human_seqcol_digests:default
	- donaldcampbelljr/human_seq_col_results:default
- Ensembl primary assembly patches (90-113, 9 total) 
	- donaldcampbelljr/ensembl_38_seqcol_digests:default
	- donaldcampbelljr/ensembl_hg38_seqcol_results:default
- NCBI GRCh genomic patches (p0-p14, 9 total), GCF  
	- donaldcampbelljr/ncbi_38_seqcol_digests:default
	- donaldcampbelljr/ncbi_38_seqcol_results:default
- NCBI p14 ->  intra patch “flavors”  
	- donaldcampbelljr/ncbi_hg38_p14_subset_digests:default
- UCSC hg19 p13 -> intra patch “flavors”
	- donaldcampbelljr/ucsc_hg19_subset_digests:default
	- see the main results pep for results

### Mus musculus

- all sampled genomes for Mus musculus (~32 genomes)
    - donaldcampbelljr/mouse_seqcol_digests:default
    - donaldcampbelljr/mouse_seq_col_results:default