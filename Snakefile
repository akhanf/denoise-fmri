from snakebids import bids

configfile: 'config.yml'

suffixes = list(config['in_images'].keys())
zip_subjects = config['zip_subjects']
zip_sessions = config['zip_sessions']

sessions = config['sessions']
if sessions == None:
    subj_wildcards = { 'subject': '{subject}'}
else:
    subj_wildcards = { 'subject': '{subject}', 'session': '{session}' }
    

rule all: 
    input: 
        expand(bids(root=config['output_dir'],desc='preproc',suffix=f'{suffix}.nii.gz',**subj_wildcards), zip, subject=zip_subjects[suffix], session=zip_sessions[suffix]) for suffix in suffixes


def get_input_images(wildcards):
    if 'session' in wildcards.keys():
        return [entry['path'] for entry in config['in_images'][wildcards.suffix][wildcards.subject][wildcards.session]]
    else:
        return [entry['path'] for entry in config['in_images'][wildcards.suffix][wildcards.subject]]


rule import_image:
    input: get_input_images
    output: bids(root=config['output_dir'],desc='preproc',suffix='{suffix}.nii.gz', **subj_wildcards)
