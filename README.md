# refgenomes.databio.org server overview

This repository contains the files to build and archive genome assets to serve with [refgenieserver](https://github.com/refgenie/refgenieserver) at http://refgenomes.databio.org. 

The whole process is scripted, starting from this repository. From here, we do this basic workflow:

1. Download raw input files for assets (FASTA files, GTF files etc.)
2. Build assets with `refgenie build` in a local refgenie instance
3. Archive assets with `refgenieserver archive`
4. Deploy assets to active server on AWS.


# Adding an asset to this server

## Overview of metadata structure

The metadata is located in the [asset_pep](asset_pep) folder, which contains a [PEP](https://pep.databio.org) with metadata for each asset. The contents are:

- `assets.csv` - The primary sample_table. Each each row is an asset. 
- `recipe_inputs.csv` - The subsample_table. This provides a way to define each individual value passed to any of the 3 arguments of the `refgenie build` command: `--assets`, `--params`, and `--files`. 
- `refgenie_build_cfg.yaml` -- config file that defines a subproject (which is used to download the input data) and additional project settings.

## Step 1: Add the asset to the asset table.

To add an asset, you will need to add a row in `assets.csv`. Follow these directions:

- `genome` - the human-readable genome (namespace) you want to serve this asset under
- `asset` - the human-readble asset name you want to serve this asset under. It is identical to the asset recipe. Use `refgenie list` to see [available recipes](http://refgenie.databio.org/en/latest/build/)

Your asset will be retrievable from the server with `refgenie pull {genome}/{asset_name}`.

## Step 2: Add any required inputs to the recipe_inputs table

Next, we need to add the source for each item required by your recipe. You can see what the recipe requires by using `-q` or `--requirements`, like this: `refgenie build {genome}/{recipe} -q`. If your recipe doesn't require any inputs, then you're done. If it requires any inputs (which can be one or more of the following: *assets*, *files*, *parameters*), then you need to specify these in the `recipe_inputs.csv` table.

For each required input, you add a row to `recipe_inputs.csv`. Follow these directions:
- `sample_name` - must match the `genome` and `asset` value in the `assets.csv` file. Format it this way: `<genome>-<asset>`. This is how we match inputs to assets.

Next you will need to fill in 3 columns:
- `input_type` which is one of the following: *files*, *params* or *assets*
- `intput_id` must match the recipe requirement. Again, use `refgenie build <genome>/<asset> -q` to learn the ids
- `input_value` value for the input, e.g. URL in case of *files*

## Step 3: See if you did it well!

**Validate the PEP with [`eido`](http://eido.databio.org/en/latest/)**

The command below validates the PEP aginst a remote schema. Any PEP issues will result in a `ValidationError`:

```
eido validate refgenie_build_cfg.yaml -s http://schema.databio.org/refgenie/refgenie_build.yaml
```

# Deploying assets onto the server

## Setup

```
#export BASEDIR=$HOME/code/sandbox/refgenie_deploy
#export REFGENIE_RAW=$BASEDIR/refgenie_raw
export BASEDIR=$PROJECT/deploy/rg.databio.org
cd $BASEDIR
git clone git@github.com:refgenie/refgenomes.databio.org.git
```

GENOMES points to pipeline output (referenced in the project config)

```
export GENOMES=$BASEDIR/genomes
```

## Step 1: Download input files

Many of the assets require some input files, and we have to make sure we have those files locally. In the `recipe_inputs.csv` file, we have entered these files as remote URLs, so the first step is to download them. We have created a subproject called `getfiles` for this: To programmatically download all the files required by `refgenie build`, run from this directory using [looper](http://looper.databio.org):

```
cd refgenomes.databio.org
export REFGENIE_RAW=/project/shefflab/www/refgenie_raw
mkdir -p $REFGENIE_RAW
looper run asset_pep/refgenie_build_cfg.yaml -p local --amend getfiles
```

## Step 2: Build assets

Once files are present locally, we can run `refgenie build` on each asset specified in the sample_table (`assets.csv`):

```
export REFGENIE=$BASEDIR/refgenomes.databio.org/config/refgenie_config.yaml
looper run asset_pep/refgenie_build_cfg.yaml -p bulker_slurm --limit 5
```

This will create one job for each *asset*. Monitor job progress with: 

```
looper check asset_pep/refgenie_build_cfg.yaml  # this doesn't work?
```

## Step 3. Archive assets

Assets are built locally now, but to serve them, we must archive them using `refgenieserver`. The command is simple:

```
refgenieserver archive -c <path/to/genomes.yaml>
```

Since the archive process is generally lengthy, it makes sense to submit this job to the cluster. We can use looper to that.

```
ba
export REFGENIE_ARCHIVE=$GENOMES/archive
looper run asset_pep/refgenie_build_cfg.yaml

```

Now we'll sync to aws. Use the refgenie credentials (here added with `--profile refgenie`, which should be preconfigured with `aws configure`)

```
aws s3 sync $REFGENIE_ARCHIVE s3://cloud.databio.org/refgenie --profile refgenie
```

## Step 4. Deploy server 

```
refgenieserver serve genomes.yaml
```

Changes to the refgenie config files with automatically trigger deploy jobs to push the updates to AWS ECS. There's a workflow for each the *master* and *staging* config file; if you change one, it will automatically deploy the correct thing.

```
ga -A; gcm "Deploy to ECS"; gpoh
```
