from bids import BIDSLayout
import os
import json
import re

# use pybids to get the paths to input images, then create a input path and wildcards for each suffix type

def read_bids_tags(bids_json=None):
    if bids_json == None:
        bids_json = os.path.join(os.path.dirname(os.path.realpath(__file__)),'bids_tags.json')
    with open(bids_json, 'r') as infile:
        bids_tags = json.load(infile)
    return bids_tags

def get_input_config_from_bids(bids_layout, inputs_dict, **filters ):
    """ returns: dict with input_path and input_wildcards"""

    bids_tags = read_bids_tags()

    config = dict({'input_path': {}, 'input_zip_lists': {}, 'input_lists': {}, 'input_wildcards': {}})





    for input_name in inputs_dict.keys():
        
        imgs, = [bids_layout.get(**inputs_dict[input_name]['filters'], **filters)]
        paths = set()
        zip_lists = {}
        input_lists = {}
        wildcards = {}
        for img in imgs:
            path = img.path
            for wildcard_name in inputs_dict[input_name]['wildcards']:

                if wildcard_name in bids_tags:
                    tag = bids_tags[wildcard_name]
                else:
                    tag = wildcard_name  #if it's not in the bids_tags dictionary, then just use the name itself as the tag

                   


                #this changes e.g. sub-001 to sub-{subject} in the path (so snakemake can use the wildcards)
                if wildcard_name in img.get_entities():
                    
                    ##HACK FIX FOR acq vs acquisition etc  -- should eventually update the bids() function to also use bids_tags.json, where e.g. acquisition -> acq is defined.. -- then, can use wildcard_name instead of out_name.. 
                    if wildcard_name not in ['subject', 'session']:
                        out_name = tag
                    else:
                        out_name = wildcard_name
 
                    if out_name not in zip_lists:
                        zip_lists[out_name] = []
                        input_lists[out_name] = set()
                        wildcards[out_name] = {}
                    pattern = '{tag}-([a-zA-Z0-9]+)'.format(tag=bids_tags[wildcard_name])
                    replace = '{tag}-{{{replace}}}'.format(tag=bids_tags[wildcard_name],replace=out_name)
                    match = re.search(pattern,path)
                    replaced = re.sub(pattern,  replace , path)
                    #update the path with the {wildcards} -- uses the value from the string (not from the pybids entities), since that has issues with integer formatting (e.g. for run=01)
                    path = replaced
                    zip_lists[out_name].append(match[1])
                    input_lists[out_name].add(match[1])
                    wildcards[out_name] = f'{{{out_name}}}'
                
            paths.add(path)
        

        #now, check to see if unique
        if len(paths) > 1:
            print(f'ERROR: more than one snakemake filename for {input_name}:')
            print(f'  Either add new bids entities to {input_name} -> wildcards, or filters to narrow the search')
            print(paths)
            return None

        in_path = list(paths)[0]

        #convert sets to lists
        for key,val in input_lists.items():
            input_lists[key] = list(val)
            
            
        config['input_path'][input_name] = in_path
        config['input_zip_lists'][input_name] = zip_lists
        config['input_lists'][input_name] = input_lists
        config['input_wildcards'][input_name] = wildcards


    return config




