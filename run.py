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
from snakebids.inputs import get_input_config_from_bids

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

#create inputs config file 
# (uses pybids to search/grab, stores lookups in a config yaml for snakemake)
inputs_config = os.path.join(args.output_dir,'code',f'inputs_config.yml')



layout = BIDSLayout(args.bids_dir,derivatives=False,validate=False,index_metadata=False, **search_terms)


inputs_config_dict = get_input_config_from_bids(bids_layout=layout, inputs_dict=config['pybids_inputs'], **search_terms)

inputs_config_dict['subjects'] = layout.get_subjects(**search_terms)
inputs_config_dict['sessions'] = layout.get_sessions(**search_terms)
if len(inputs_config_dict['sessions'])  == 0:
    inputs_config_dict['subj_wildcards'] = { 'subject': '{subject}'}
else:
    inputs_config_dict['subj_wildcards'] = { 'subject': '{subject}', 'session': '{session}' }

#write updated config file
os.makedirs(os.path.dirname(inputs_config),exist_ok=True)

with open(inputs_config, 'w') as outfile:
    yaml.dump(inputs_config_dict, outfile, default_flow_style=False)




# running participant level
if args.analysis_level == "participant":

    #runs snakemake, using the workflow config and inputs config to override 
    snakefile = os.path.join(snakemake_dir,'workflow/Snakefile')
    
    if args.use_snakemake_api:
        snakemake(snakefile,configfiles=[workflow_config, inputs_config], workdir=args.output_dir, dryrun=True, debug_dag=True)
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

