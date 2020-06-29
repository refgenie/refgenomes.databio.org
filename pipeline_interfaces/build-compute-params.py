#!/usr/bin/env python3

import json
from argparse import ArgumentParser

parser = ArgumentParser(description="Refgenie build params")

parser.add_argument("-s", "--size", help="size", required=False)
parser.add_argument("-a", "--asset", type=str, help="asset", required=True)

args = parser.parse_args()

params =   {
    "bulker_crate": "databio/refgenie:0.7.6",
    "mem": "24000",
    "cores": "1",
    "partition": "largemem",
    "time": "04:00:00"}

# These are ones that go quick
fast_assets = ["fasta", "gencode_gtf", "ensembl_gtf", "ensembl_rb",
               "refgene_anno", "dbnsfp", "fasta_txome"]

slow_assets = ["suffixerator_index", "salmon_partial_sa_index"]


if args.asset in fast_assets:
    params['time'] = "00:30:00"

if args.asset in slow_assets:
    params['time'] = "8:00:00"

if args.asset == 'suffixerator_index':
    params['mem'] = "32000"

if args.asset == 'bowtie2_index':
    params['mem'] = "64000"

if args.asset == 'bismark_bt2_index':
    params['mem'] = "64000"
    params["time"] = "08:00:00"

if args.asset == 'bismark_bt1_index':
    params['mem'] = "64000"
    params["time"] = "08:00:00"

if args.asset == 'salmon_partial_sa_index':
    params['mem'] = "96000"

if args.asset == 'salmon_sa_index':
    params['mem'] = "72000"

if args.asset == 'star_index':
    params['mem'] = "64000"


y = json.dumps(params)

print(y)