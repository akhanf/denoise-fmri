from nilearn.image import clean_img 
import pandas as pd
import numpy as np
from snakebids import write_derivative_json

#print('params')
#print(snakemake.params)

#print('confounds to use:')
#print(snakemake.params.confounds_to_use)

#read confounds table from fmriprep into a pandas dataframe
df = pd.read_table(snakemake.input.confounds_tsv)

#get specified confounds from the dataframe
confounds = df[snakemake.params.confounds_to_use].to_numpy()

#print('before removing non-finite')
#print(confounds)
#replace non-finite values with zeros in the confounds
confounds[np.logical_not(np.isfinite(confounds))] = 0


#use nilearn to clean
cleaned = clean_img(snakemake.input.nii, detrend=snakemake.params.detrend, standardize=snakemake.params.standardize, confounds=confounds, mask_img=snakemake.input.mask_nii)

#save to nifti
cleaned.to_filename(snakemake.output.nii)

#pass along the derivatives sidecar, with addition of Description and Sources 
write_derivative_json(snakemake, Description=f'Denoising using confound grouping {snakemake.params.confounds_name}')


