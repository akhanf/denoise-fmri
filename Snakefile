from snakebids import bids

configfile: 'config.yml'

#get all possible suffixes
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


def get_input_image(wildcards):
    if 'session' in wildcards.keys():
        return config['in_images'][wildcards.suffix][wildcards.subject][wildcards.session]['path']
    else:
        return config['in_images'][wildcards.suffix][wildcards.subject]['path']
    

rule import_image:
    input: get_input_image
    output: bids(root=config['output_dir'],desc='preproc',suffix='{suffix}.nii.gz', **subj_wildcards)
