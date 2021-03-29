#!/usr/bin/env python

import sys
import tqdm
import requests
import pandas as pd

with open(sys.argv[1]) as fh:
    run_accessions = [x.strip() for x in fh]

# break into 500 accession chunks to more efficiently query the API
api_responses = []

for accession in tqdm.tqdm(run_accessions):
    api_response = requests.get(f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_md5")
    print(api_response.text.split('\n')[1])

