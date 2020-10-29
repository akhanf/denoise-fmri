from nilearn.image import clean_img 
import pandas as pd

#read confounds table from fmriprep into a pandas dataframe
df = pd.read_table(snakemake.input.confounds_tsv)

#get specified confounds from the dataframe
confounds = df[snakemake.params.confounds_to_use].to_numpy()

#use nilearn to clean
cleaned = clean_img(snakemake.input.nii, detrend=True, standardize=True, confounds=confounds, mask_img=snakemake.input.mask_nii)

#save to nifti
cleaned.to_filename(snakemake.output.denoised)
