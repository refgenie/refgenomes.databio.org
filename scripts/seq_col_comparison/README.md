# Sequence Collection Comparison Pipeline

### Sample Level Pipeline

For each sample in your PEP:
- create sequence_collection dictionary and dump it to a json
- create chrom.sizes files from the sequence_collection dict

### Project Level Pipeline

Calculate stats for all samples in the PEP including, jaccard similarities for names, lengths, name_len_pairs, sequences.

Future optimization: this project-level pipeline is currently inefficient and take several minutes for N>45

### Plotting statistics figures
After running project pipeline and calculating statistics, you can plot figures. 

You need to provide the PEPhub path for the PEP with stats results, species (for figure titles), PEPhub path that links sample_name with authority, path to local directory with jsons of sequence collection information.

Example command
```
python3 stats_graphs.py /home/drc/Downloads/refgenomes_pics_test/22may2025/mus_musculus/ donaldcampbelljr/mouse_seq_col_results:default mus_musculus /home/drc/Downloads/mouse_jsons_from_rivanna/json/ donaldcampbelljr/mouse_seqcol_digests:default 
```

There is also an R script, `test_plotting.R`. You must add the path to the output folder for the above python script (it produces csvs for plotting)

#### plotting clustered sequence frequency
- run script `cluster_seqcounts.py` and then point `cluster_sequence_plots.R` to its output and run to produce a figure of sequences clustered by frequency.