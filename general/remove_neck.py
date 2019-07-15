#!/usr/bin/python

import os
import nibabel as nib
import numpy as np

# Load in structural
filename = '/scratch/groups/Projects/P1323/fMRI/R4127/Alignment/High_Res/nu.nii.gz'
dirname, basename = os.path.split(filename)
nii = nib.load(filename)

# Get various useful bits
data = nii.get_data() # actual contents of volume
xdim, ydim, zdim = data.shape
header = nii.get_header() # useful later
affine = nii.get_affine() # useful later

# Apply mask
cut_data = data.copy()

#cut_data[:x:y:z] = 0

cut_data[:,:,:74] = 0

# Save out
nib.Nifti1Image(cut_data, header = header, affine = affine).to_filename(os.path.join(dirname, 'Fixed_' + basename))

# TO RUN THIS: Go into terminal, cd into the 'toolbox' dir that this code is saved in. Type ipython then type run and drag and drop the script.
