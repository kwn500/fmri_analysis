#!/bin/bash

## run bet on nu from freesurfer - non-skull stripped.
#bet nu.nii.gz bet_nu -R
#??? skull stripping???


# Bias correct the flair using a hacked version of FAST (actually fs_anat)
/scratch/groups/Projects/P1323/Code/fsl_prep/fsl_anat_local -i '/scratch/groups/Projects/P1323/fMRI/@RNUMBER@/Alignment/Inplane/@INPLANEANATFILE@' -o /scratch/groups/Projects/P1323/fMRI/@RNUMBER@/Alignment/corrFlair -t T1 --clobber --noreorient --nocrop --noreg --nononlinreg --noseg --nosubcortseg --nobet --nosearch 

# try to align the bias-corrected non-skull stripped T1 flair to the equivalent whole-brain T1
flirt -in /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/corrFlair.anat/T1_biascorr.nii.gz -ref /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/High_Res/nu.nii.gz -out /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/initial_highres2highres -omat /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/initial_highres2highres.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

# Bias correct the example_func using a hacked version of FAST (actually fs_anat)
/scratch/groups/Projects/P1323/Code/fsl_prep/fsl_anat_local -i '/scratch/groups/Projects/P1323/fMRI/FeatOutput/@CONDITION@_Analysis/Level_1/FEAT_Final_@RNUMBER@_Run@RUNNUMBER@.feat/example_func.nii.gz' -o /scratch/groups/Projects/P1323/fMRI/@RNUMBER@/Alignment/corrExampleFunc -t T2 --clobber --noreorient --nocrop --noreg --nononlinreg --noseg --nosubcortseg --nobet --nosearch 

# try to align the bias-corrected non-skull stripped example func to the equivalent flair
flirt -in /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/corrExampleFunc.anat/T2_biascorr.nii.gz -ref /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/corrFlair.anat/T1_biascorr.nii.gz -out /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/example_func2initial_highres.nii.gz -omat /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/example_func2initial_highres.mat -bins 256 -cost mutualinfo -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

#and align the highres to standard
flirt -in /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/High_Res/@HIGHRESANATFILE@ -ref /usr/share/fsl-5.0/data/standard/MNI152_T1_2mm -out /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/highres2standard -omat /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/highres2standard.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear

#copy the 3 transformation matrices to  the reg directory of the feat session
cp /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/example_func2initial_highres.mat /scratch/groups/Projects/P1323/fMRI/FeatOutput/@CONDITION@_Analysis/Level_1/FEAT_Final_@RNUMBER@_Run@RUNNUMBER@.feat/reg/
cp /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/initial_highres2highres.mat /scratch/groups/Projects/P1323/fMRI/FeatOutput/@CONDITION@_Analysis/Level_1/FEAT_Final_@RNUMBER@_Run@RUNNUMBER@.feat/reg/
cp /scratch/sg3/P1323/fMRI/@RNUMBER@/Alignment/highres2standard.mat /scratch/groups/Projects/P1323/fMRI/FeatOutput/@CONDITION@_Analysis/Level_1/FEAT_Final_@RNUMBER@_Run@RUNNUMBER@.feat/reg/

updatefeatreg /scratch/groups/Projects/P1323/fMRI/FeatOutput/@CONDITION@_Analysis/Level_1/FEAT_Final_@RNUMBER@_Run@RUNNUMBER@.feat
featregapply /scratch/groups/Projects/P1323/fMRI/FeatOutput/@CONDITION@_Analysis/Level_1/FEAT_Final_@RNUMBER@_Run@RUNNUMBER@.feat
