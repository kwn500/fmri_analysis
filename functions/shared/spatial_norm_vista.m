ni = readFileNifti('/scratch/groups/Projects/P1361/fMRI/R5006/anatomy/t1.nii.gz');
[sn,template,inv] = mrAnatComputeSpmSpatialNorm(double(ni.data), ni.qto_xyz, 'MNI_T1');
invLUT = rmfield(inv,{'deformX','deformY','deformZ'});
save('/scratch/groups/Projects/P1361/fMRI/R5006/anatomy/t1_sn.mat','sn','invLUT');