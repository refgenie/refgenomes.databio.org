#!/usr/bin/env python3

import json
from argparse import ArgumentParser

parser = ArgumentParser(description="Refgenie build params")

parser.add_argument("-s", "--size", help="size", required=False)
parser.add_argument("-a", "--asset", type=str, help="asset", required=True)
parser.add_argument("-g", "--genome", type=str, help="genome", required=True)

args = parser.parse_args()

params =   {
    "bulker_crate": "databio/refgenie:0.7.6",
    "mem": "24000",
    "cores": "1",
    "partition": "largemem",
    "time": "04:00:00"}

# These are ones that go quick
fast_assets = ["fasta", "gencode_gtf", "ensembl_gtf", "ensembl_rb", "feat_annotation",
               "refgene_anno", "fasta_txome"]

slow_assets = ["bismark_bt2_index", "bismark_bt1_index", "salmon_partial_sa_index"]


if args.asset in fast_assets:
    params['time'] = "00:30:00"
    params['partition'] = "standard"
    params['mem'] = "6000"

if args.asset in slow_assets:
    params['time'] = "8:00:00"

if args.asset == 'suffixerator_index':
    params['mem'] = "32000"

if args.asset == 'bowtie2_index':
    params['mem'] = "64000"

if args.asset == 'bismark_bt2_index':
    params['mem'] = "64000"

if args.asset == 'bismark_bt1_index':
    params['mem'] = "64000"

if args.asset == 'salmon_partial_sa_index':
    params['mem'] = "112000"
    params['time'] = "6:00:00"
    params['cores'] = "8"

if args.asset == 'dbsfnp':
    params['time'] = "12:00:00"

if args.asset == 'salmon_sa_index':
    params['mem'] = "72000"

if args.asset == 'star_index':
    params['mem'] = "64000"


y = json.dumps({"compute": params})

print(y)
