#!/usr/bin/env python3
import os
import subprocess
import sys
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

    import yaml
#    from datetime import datetime
#    date_string = datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')

#    config_file = os.path.join(args.output_dir,'code',f'config_{date_string}.yml')
#    os.makedirs(os.path.dirname(config_file),exist_ok=True)

    config_file = 'inputs_config.yml'

    with open(config_file, 'w') as outfile:
        yaml.dump(config_dict, outfile, default_flow_style=False)

    return config_file



#start of main

parser =  get_parser()
args = parser.parse_args()

if not args.skip_bids_validator:
    run('bids-validator %s'%args.bids_dir)



#default search term is by suffix
search_terms = {'suffix': args.suffix}

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
# for all subjects
else:
    subjects = layout.get_subjects(**search_terms)
    sessions = layout.get_sessions(**search_terms)

if len(sessions) == 0:
    sessions = None



# what to pass to snakemake?
#  set --config variables?  could have a schema defined as a standard.. 
#   e.g.:  subjects (list), input_files (dict, indexed by subject), output_folder (string), 
#
# with API, can either pass config dict (but must be flat! ie cannot be nested), or make a config json/yaml, and pass that as a config file.. 
#   latter is probably better, as there is provenance too with the configfile.. 


#different ways reading input files for a given subject:
# 1) take list of images (from particular suffix), to be merged in the pipeline  (iterate over subject)        ** this will be similar to 2)
#       in this case, want the input image list to be indexed by subject
#
# 2) take a single image for each subject/session (iterate over subject and session)           ** implement this for now
#       in this case, want the input image to be indexed by subject (and optional session) 
#
# 3) iterate over all found images (sessions, acquisitions etc)
#       in this case, want tthe input image to be indexed by all the entities (to make it unique)       ** this one will take thought.. may not be practically used very often.. 
#
# hash entities into a bids string -> this is unique -- this can be the wildcard??
# 1-1 mapping between input files and bids string hashes
# e.g.:   sub-01/sub-01_run-01.nii.gz,  sub-01/sub-01_run-02.nii.gz
# wildcard could be the entire bids-string, and could have a get_entities() function that pulls out the key/vals

    

# running participant level
if args.analysis_level == "participant":


    # use pybids to get the paths to input images, then create a dict for each suffix type
    #   which is indexed by subject, and contains the path and entities
    
    config_dict = {'subjects': subjects, 'sessions': sessions,'output_dir': args.output_dir, 'zip_subjects':{}, 'zip_sessions':{}, 'in_images': {}}
    for suffix in args.suffix:
        print(suffix)
        config_dict['zip_subjects'][suffix] = []
        config_dict['zip_sessions'][suffix] = []
        config_dict['in_images'][suffix] = {}

        #refine to single search term
        search_terms['suffix'] = suffix

        
        #for this type of image, want no more than 1 per session    
        all_imgs = []
        for subject in subjects:
            if sessions == None:
                imgs = layout.get(subject=subject,extension='nii.gz',**search_terms)
                if len(imgs) > 1:
                    print('ERROR: more than one image matched for sub-{subject} please specify additional search terms to avoid ambiguity')
                    sys.exit()
                elif len(imgs) == 1:
                    all_imgs.append(*imgs)

            else:
                for session in sessions:
                    imgs = layout.get(subject=subject,session=session,extension='nii.gz',**search_terms)
                    if len(imgs) > 1:
                        print('ERROR: more than one image matched for sub-{subject}/ses-{session}, please specify additional search terms to avoid ambiguity')
                        sys.exit()
                    elif len(imgs) == 1:
                        all_imgs.append(*imgs)

        #iterate through images, adding to the config file that snakemake will use
        for img in all_imgs:

            entities = img.get_entities()
            subject = entities['subject']
            suffix = entities['suffix']

            if 'session' in entities.keys(): 
                session = entities['session']
                if subject in config_dict['in_images'][suffix]: 
                    config_dict['in_images'][suffix][subject].update(  {session: {'path': img.path, 'entities': {**entities} }} )
                else: 
                    config_dict['in_images'][suffix][subject] = {session: {'path': img.path, 'entities': {**entities} }}
                config_dict['zip_subjects'][suffix].append(subject)
                config_dict['zip_sessions'][suffix].append(session)

            else:
                if sessions == None:
                    config_dict['in_images'][suffix][subject] = {'path': img.path, 'entities': {**entities} }
                else:
                    print(f'ERROR: {img.path}, sessions must be defined for all images if sessions are used')
                    sys.exit()
                config_dict['zip_subjects'][suffix].append(subject)
            

    # create a config file, and pass that:
    config_file = create_yaml_cfg(config_dict)

    snakemake('Snakefile',configfiles=[config_file], dryrun=True, printshellcmds=True)

    

# running group level
elif args.analysis_level == "group":
    print('insert report generation here!')

