#!/usr/bin/env python3

from snakebids.app import SnakeBidsApp

#app = SnakeBidsApp()


app = SnakeBidsApp(snakebids_config='/scratch/akhanf/bids_app_smk/snakebids-app/config/snakebids.yml',
                    snakefile='/scratch/akhanf/bids_app_smk/snakebids-app/workflow/Snakefile')
app.run_snakemake()
