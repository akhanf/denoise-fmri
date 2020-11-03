#!/usr/bin/env python3

from snakebids.app import SnakeBidsApp


app = SnakeBidsApp(snakebids_config='./config/snakebids.yml',
                    snakefile='./workflow/Snakefile')
app.run_snakemake()
