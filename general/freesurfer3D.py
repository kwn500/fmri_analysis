#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 13:22:37 2019

@author: kwn500
"""

from surfer import Brain 
import os

subject_id = 'fsaverage'
hemi = 'lh'
surf = 'inflated'
cope_filepath = '/scratch/groups/Projects/P1361/fMRI/fsl_output/level3/freesurfer/'
cope_output = cope_filepath+'output/'
cope_file = 'shape_versus_passive_zstat_lh'

os.environ['SUBJECTS_DIR'] = "/usr/share/freesurfer-data/subjects"
brain = Brain(subject_id, hemi, surf)
brain.add_overlay(cope_filepath+cope_file+'.mgh', 2.3, 8)


brain.save_imageset(cope_output+cope_file, ['medial', 'lateral', 'rostral', 'caudal', 'dorsal', 'ventral'], 'png', [1,2,3,4,5,6])
brain.close()