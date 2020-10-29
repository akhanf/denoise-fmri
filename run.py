#!/usr/bin/env python3
import os
import subprocess
import sys
import yaml
from bids import BIDSLayout
import bids
from glob import glob
from snakemake import snakemake
import argparse
import json

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

#get path of this file
snakemake_dir = os.path.dirname(os.path.realpath(__file__))

#replace with proper name of pipeline here
parser = argparse.ArgumentParser(description='snakebids-app')

#load up workflow config file
snakebids_config = os.path.join(snakemake_dir,'cfg','snakebids.yml')
with open(snakebids_config, 'r') as infile:
    config = yaml.load(infile, Loader=yaml.FullLoader)

#update the parser with config options
for name, parse_args in config['parse_args'].items():
    parser.add_argument(name, **parse_args)

all_args = parser.parse_known_args()

args = all_args[0]
snakemake_args = all_args[1]

#add arguments to config
config.update(args.__dict__)

search_terms = dict()

#add optional search terms
if args.participant_label is not None:
    search_terms['subject'] = args.participant_label
if args.session is not None:
    search_terms['session'] = args.session

# still  need to come up with ways to specify modality-specific search terms from command-line (other than --config ...)
"""
if args.acq is not None:
    search_terms['acquisition'] = args.acq
if args.run is not None:
    search_terms['run'] = args.run
"""

#override bids_dir, output_dir, search_terms in snakebids config
config['search_terms'] = search_terms
config['bids_dir'] = os.path.realpath(args.bids_dir)
config['output_dir'] = os.path.realpath(args.output_dir)


#create an updated snakebids config file
updated_config = os.path.join(config['output_dir'],'cfg','snakebids.yml')

#write it to file
os.makedirs(os.path.dirname(updated_config),exist_ok=True)
with open(updated_config, 'w') as outfile:
    yaml.dump(config, outfile, default_flow_style=False)

# run snakemake with that config


#workflow snakefile will read snakebids config, create inputs_config, read that in


# running participant level
if args.analysis_level == "participant":

    #runs snakemake, using the workflow config and inputs config to override 
    snakefile = os.path.join(snakemake_dir,'workflow/Snakefile')
    
    if args.use_snakemake_api:
        snakemake(snakefile,configfiles=[updated_config], workdir=args.output_dir, dryrun=True, debug_dag=True)
    else:
        #run snakemake command-line (passing any leftover args from argparse)
        snakemake_cmd_list = ['snakemake',
                                f'--snakefile {snakefile}',
                                f"--directory {config['output_dir']}",
                                f'--configfile {updated_config} ',
                                *snakemake_args]

        snakemake_cmd = ' '.join(snakemake_cmd_list)
        run(snakemake_cmd)


# running group level
elif args.analysis_level == "group":
    print('insert report generation here!')

