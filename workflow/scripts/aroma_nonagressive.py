from nilearn.image import smooth_img 
from snakebids import write_derivative_json

#use nilearn to smooth
smoothed = smooth_img(snakemake.input.nii, fwhm=snakemake.params.fwhm)

#save to niftii
smoothed.to_filename(snakemake.output.nii)

#pass along the derivatives sidecar, with addition of Description and Sources 
write_derivative_json(snakemake, Description=f'Smoothing with FWHM={snakemake.params.fwhm}')

"""
#save sidecar json
import json
with open(snakemake.input.json,'r') as f:
  sidecar = json.load(f)

sidecar.update({'Description': f'Smoothing with FWHM={snakemake.params.fwhm}',
                'Sources': [snakemake.input] })

with open(snakemake.output.json, 'w') as outfile:
    json.dump(sidecar, outfile,indent=4)

"""
