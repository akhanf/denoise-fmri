from bids import BIDSLayout
import yaml
import os


# use pybids to get the paths to input images, then create a dict for each suffix type
#   which is indexed by subject, and contains the path and entities

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



def create_bids_input_config(bids_dir, suffixes, out_config_yml, participant_label=None, search_terms={}):

    #start by getting bids layout from pybids
    layout = BIDSLayout(bids_dir,derivatives=False,validate=False,index_metadata=False)

    # only for a subset of subjects
    if participant_label:
        subjects = participant_label
        sessions = layout.get_sessions(subject=subjects,**search_terms)
    # get all subjects
    else:
        subjects = layout.get_subjects(**search_terms)
        sessions = layout.get_sessions(**search_terms)

    if len(sessions) == 0:
        sessions = None

    # initialize the config_dict 
    config_dict = {'subjects': subjects, 'sessions': sessions, 
                    'zip_subjects':{}, 'zip_sessions':{}, 'in_images': {}, 'search_terms': {}}


    #populate config file with the input images
    for suffix in suffixes:
        config_dict = add_images_to_config(config_dict=config_dict, suffix=suffix, search_terms=search_terms, layout=layout)

#    config_file = os.path.join(args.output_dir,'code',f'inputs_config.yml')

    os.makedirs(os.path.dirname(out_config_yml),exist_ok=True)

    with open(out_config_yml, 'w') as outfile:
        yaml.dump(config_dict, outfile, default_flow_style=False)



