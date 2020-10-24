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



def create_yaml_cfg( config_dict):
    """Create config yml using a dict, and return path to it in the output folder"""


    config_file = os.path.join(args.output_dir,'code',f'inputs_config.yml')
    os.makedirs(os.path.dirname(config_file),exist_ok=True)

    with open(config_file, 'w') as outfile:
        yaml.dump(config_dict, outfile, default_flow_style=False)

    return config_file




def add_images_to_config(config_dict, suffix, search_terms, layout):
    """ returns: config_dict: updated config with fields updated
                -zip_subjects, zip_sessions, in_images, search_terms  """
    
    subjects = config_dict['subjects']
    sessions = config_dict['sessions']

    #add suffix to search terms
    search_terms['suffix'] = suffix

    #initialize the new entries to config_dict (will be updated in loop below)
    config_dict['zip_subjects'][suffix] = []
    config_dict['zip_sessions'][suffix] = []
    config_dict['in_images'][suffix] = {}
    config_dict['search_terms'][suffix] = {**search_terms}
    
    all_imgs = []

    #first, get the list of bidsimages
    for subject in subjects:
        if sessions == None:
            imgs = layout.get(subject=subject,extension='nii.gz',**search_terms)
            if imgs is not None:
                all_imgs = all_imgs + imgs
                config_dict['in_images'][suffix][subject] = []
                config_dict['zip_subjects'][suffix].append(subject)
        else:
            for session in sessions:
                imgs = layout.get(subject=subject,session=session,extension='nii.gz',**search_terms)
                if imgs is not None:
                    all_imgs = all_imgs + imgs
                    config_dict['in_images'][suffix][subject] = {session: [] }
                    config_dict['zip_subjects'][suffix].append(subject)
                    config_dict['zip_sessions'][suffix].append(session)


    #iterate through images, adding to the config file that snakemake will use
    for img in all_imgs:

        entities = img.get_entities()
        subject = entities['subject']
        suffix = entities['suffix']

        if 'session' in entities.keys(): 
            session = entities['session']
            config_dict['in_images'][suffix][subject][session].append( {'path': img.path, 'entities': {**entities} } )
        else:
            if sessions == None:
                config_dict['in_images'][suffix][subject].append( {'path': img.path, 'entities': {**entities}})
            else:
                print(f'ERROR: {img.path}, sessions must be defined for all images if sessions are used')
                sys.exit()

    return config_dict





#start of main

parser =  get_parser()
args = parser.parse_args()

if not args.skip_bids_validator:
    run('bids-validator %s'%args.bids_dir)



#default search term is by suffix
search_terms = {'suffix': ['T2w', 'dwi']}

#add optional search terms
if args.session is not None:
    search_terms['session'] = args.session
if args.acq is not None:
    search_terms['acquisition'] = args.acq
if args.run is not None:
    search_terms['run'] = args.run


#start by getting bids layout from pybids
layout = BIDSLayout(args.bids_dir,derivatives=False,validate=False,index_metadata=False)

# only for a subset of subjects
if args.participant_label:
    subjects = args.participant_label
    sessions = layout.get_sessions(subject=subjects,**search_terms)
# get all subjects
else:
    subjects = layout.get_subjects(**search_terms)
    sessions = layout.get_sessions(**search_terms)

if len(sessions) == 0:
    sessions = None


    
# running participant level
if args.analysis_level == "participant":


    # use pybids to get the paths to input images, then create a dict for each suffix type
    #   which is indexed by subject, and contains the path and entities
    
    # initialize the config_dict 
    config_dict = {'subjects': subjects, 'sessions': sessions,'output_dir': args.output_dir, 'zip_subjects':{}, 'zip_sessions':{}, 'in_images': {}, 'search_terms': {}}


    #populate config file
    config_dict = add_images_to_config(config_dict=config_dict, suffix='T1w', search_terms=search_terms, layout=layout)
    config_dict = add_images_to_config(config_dict=config_dict, suffix='T2w', search_terms=search_terms, layout=layout)
    config_dict = add_images_to_config(config_dict=config_dict, suffix='dwi', search_terms=search_terms, layout=layout)


    # create yaml config file from the dict
    config_file = create_yaml_cfg(config_dict)

    # pass that file to snakemake, which supplements/overrides the workflow config file
    snakemake('Snakefile',configfiles=[config_file], dryrun=True, printshellcmds=True)

    
# running group level
elif args.analysis_level == "group":
    print('insert report generation here!')

