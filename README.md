# rg.databio.org server overview

This repository contains the files to build and archive genome assets to serve with [refgenieserver](https://github.com/refgenie/refgenieserver) at http://refgenomes.databio.org. 

The whole process is scripted, starting from this repository. From here, we do this basic workflow:

1. Download raw input files for assets (FASTA files, GTF files etc.)
2. Configure refgenie
3. Build assets with `refgenie build` in a local refgenie instance
4. Archive assets with `refgenieserver archive`
5. Upload archives to S3
6. Deploy assets to active server on AWS.


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

In this guide we'll use environment variables to keep track of where stuff goes.

- `BASEDIR` points to our parent folder where we'll do all the building/archiving
- `GENOMES` points to pipeline output (referenced in the project config)
- `REFGENIE_RAW` points to a folder where the downloaded raw files are kept
- `REFGENIE` points to the refgenie config file
- `REFGENIE_ARCHIVE` points to the location where we'll store the actual archives

```
export SERVERNAME=refgenomes.databio.org
export BASEDIR=$PROJECT/deploy/$SERVERNAME
export GENOMES=$BASEDIR/genomes
export REFGENIE_RAW=/project/shefflab/www/refgenie_$SERVERNAME
export REFGENIE=$BASEDIR/$SERVERNAME/config/refgenie_config.yaml
export REFGENIE_ARCHIVE=$GENOMES/archive
mkdir $BASEDIR
cd $BASEDIR
```

To start, clone this repository:

```
git clone git@github.com:refgenie/$SERVERNAME.git
```

## Step 1: Download input files

Many of the assets require some input files, and we have to make sure we have those files locally. In the `recipe_inputs.csv` file, we have entered these files as remote URLs, so the first step is to download them. We have created a subproject called `getfiles` for this: To programmatically download all the files required by `refgenie build`, run from this directory using [looper](http://looper.databio.org):

```
cd $SERVERNAME
mkdir -p $REFGENIE_RAW
looper run asset_pep/refgenie_build_cfg.yaml -p local --amend getfiles --sel-attr asset --sel-incl fasta
```

Check the status with `looper check --use-pipesat`:

*`--use-pipesat` option is required as of early 2021. Might not be required if you're running later on.*

```
looper check asset_pep/refgenie_build_cfg.yaml --amend getfiles --sel-attr asset --sel-incl fasta --use-pipestat
```

## Step 2: Refgenie genome configuration file initialization

This repository comes with files genome cofiguration file already defined in [`\config`](config) directory, but if you have not initialized refgenie yet or want to start over, then first you can initialize the config like this:

```
refgenie init -c config/refgenie_config.yaml -f $GENOMES -u http://awspds.refgenie.databio.org/refgenomes.databio.org/ -a $GENOMES/archive -b refgenie_config_archive.yaml
```

## Step 3: Build assets

Once files are present locally, we can run `refgenie build` on each asset specified in the sample_table (`assets.csv`). We have to submit fasta assets first:

```
looper run asset_pep/refgenie_build_cfg.yaml -p bulker_slurm --sel-attr asset --sel-incl fasta
```

This will create one job for each *asset*. Monitor job progress with: 

```
grep CANCELLED ../genomes/submission/*.log
ll ../genomes/submission/*.log
grep error ../genomes/submission/*.log
grep maximum ../genomes/submission/*.log

ll ../genomes/data/*/*/*/_refgenie_build/*.flag
ll ../genomes/data/*/*/*/_refgenie_build/*failed.flag
ll ../genomes/data/*/*/*/_refgenie_build/*completed.flag
ll ../genomes/data/*/*/*/_refgenie_build/*running.flag
ll ../genomes/data/*/*/*/_refgenie_build/*completed.flag | wc -l
cat ../genomes/submission/*.log
```

To run all the asset types:

```
looper run asset_pep/refgenie_build_cfg.yaml -p bulker_slurm
```

## Step 4. Archive assets

Assets are built locally now, but to serve them, we must archive them using `refgenieserver`. The general command is `refgenieserver archive -c <path/to/genomes.yaml>`. Since the archive process is generally lengthy, it makes sense to submit this job to a cluster. We can use looper to do that. 

To start over completely, remove the archive config file with: 

``` 
rm config/refgenie_config_archive.yaml
```

Then submit the archiving jobs with `looper run`

```
looper run asset_pep/refgenieserver_archive_cfg.yaml -p bulker_local --sel-attr asset --sel-incl fasta
```

Check progress with:

```
ll ../genomes/archive_logs/submission/*.log
grep Wait ../genomes/archive_logs/submission/*.log
grep Error ../genomes/archive_logs/submission/*.log
cat ../genomes/archive_logs/submission/*.log
```

## Step 5. Upload archives to S3

Now the archives should be built, so we'll sync them to AWS. Use the refgenie credentials (here added with `--profile refgenie`, which should be preconfigured with `aws configure`)


```
aws s3 sync $REFGENIE_ARCHIVE s3://awspds.refgenie.databio.org/refgenomes.databio.org/ --profile refgenie
```

## Step 6. Deploy server 

Now everything is ready to deploy. If using refgenieserver directly, you'll run `refgenieserver serve config/refgenieserver_archive_cfg`. We're hosting this repository on AWS and use GitHub Actions to trigger  trigger deploy jobs to push the updates to AWS ECS whenever a change is detected in the config file. 

```
ga -A; gcm "Deploy to ECS"; gpoh
```
