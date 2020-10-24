from snakebids import bids
from snakebids.inputs import create_bids_input_config

#only load config file if not already loaded by snakemake API (in run.py)
if len(config) == 0:
    #load main config file
    configfile: 'cfg/config.yml'

    #generate input config file from bids data
    create_bids_input_config(bids_dir=config['bids_dir'], suffixes=config['in_suffixes'], out_config_yml='cfg/inputs_config.yml') 
    
    #read that new config file
    configfile: 'cfg/inputs_config.yml'


suffixes = config['in_suffixes']
zip_subjects = config['zip_subjects']
zip_sessions = config['zip_sessions']
subjects = config['subjects']
sessions = config['sessions']
if sessions == None:
    subj_wildcards = { 'subject': '{subject}'}
else:
    subj_wildcards = { 'subject': '{subject}', 'session': '{session}' }
    
rule all:
    #input: expand(bids(root=config['output_dir'],desc='bet',suffix='T1w.nii.gz', **subj_wildcards), zip, subject=zip_subjects['T1w'], session=zip_sessions['T1w'])
    input: expand(bids(root=config['output_dir'],desc='bet',suffix='T1w.nii.gz', **subj_wildcards), subject=subjects,session=sessions)




rule all_imported: 
    input: 
        expand( bids(root=config['output_dir'],suffix=f'{suffix}.nii.gz',**subj_wildcards),
                zip, subject=zip_subjects[suffix], session=zip_sessions[suffix]) for suffix in suffixes


def get_input_images(wildcards):
    """returns: list of input images from the dict created by pybids"""

    if 'session' in wildcards.keys():
        return [entry['path'] for entry in config['in_images'][wildcards.suffix][wildcards.subject][wildcards.session]]
    else:
        return [entry['path'] for entry in config['in_images'][wildcards.suffix][wildcards.subject]]


rule import_image:
    input: get_input_images 
    output: bids(root=config['output_dir'],suffix='{suffix}.nii.gz', **subj_wildcards)
    shell:
        'cp {input} {output}'

rule bet_t1w:
    input: bids(root=config['output_dir'],suffix='T1w.nii.gz', **subj_wildcards)
    output: bids(root=config['output_dir'],desc='bet',suffix='T1w.nii.gz', **subj_wildcards)
    shell: 
        'bet {input} {output}'

