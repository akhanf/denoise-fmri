#!/usr/bin/env python3
import os
import subprocess
import sys
import yaml
from bids import BIDSLayout
import bids
from glob import glob
from snakemake import snakemake
from parse import get_parser
import json
from snakebids.inputs import create_bids_input_config

bids.config.set_option('extension_initial_dot', True)

def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)




#start of main

parser =  get_parser()
all_args = parser.parse_known_args()

args = all_args[0]
snakemake_args = all_args[1]


#for running snakemake
snakemake_dir = os.path.dirname(os.path.realpath(__file__))


#load up workflow config file
workflow_config = os.path.join(snakemake_dir,'cfg','config.yml')
with open(workflow_config, 'r') as infile:
    config = yaml.load(infile, Loader=yaml.FullLoader)




#default search term is by suffix
search_terms = {'suffix': ['T2w', 'dwi']}

#add optional search terms
if args.session is not None:
    search_terms['session'] = args.session
if args.acq is not None:
    search_terms['acquisition'] = args.acq
if args.run is not None:
    search_terms['run'] = args.run


#create inputs config file 
# (uses pybids to search/grab, stores lookups in a config yaml for snakemake)
inputs_config = os.path.join(args.output_dir,'code',f'inputs_config.yml')

create_bids_input_config(bids_dir=args.bids_dir, suffixes=config['in_suffixes'], 
                        out_config_yml=inputs_config, 
                        participant_label=args.participant_label,
                        search_terms=search_terms)


# running participant level
if args.analysis_level == "participant":

    #runs snakemake, using the workflow config and inputs config to override 
    snakefile = os.path.join(snakemake_dir,'Snakefile')
    
    if args.use_snakemake_api:
        snakemake(snakefile,configfiles=[workflow_config, inputs_config], workdir=args.output_dir)
    else:
        #run snakemake command-line (passing any leftover args from argparse)
        snakemake_cmd_list = ['snakemake',
                                f'--snakefile {snakefile}',
                                f'--directory {args.output_dir}',
                                f'--configfile {workflow_config} {inputs_config}',
                                *snakemake_args]

        snakemake_cmd = ' '.join(snakemake_cmd_list)
        run(snakemake_cmd)


# running group level
elif args.analysis_level == "group":
    print('insert report generation here!')

