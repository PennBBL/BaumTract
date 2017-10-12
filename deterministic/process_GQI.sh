#!/bin/sh

##################################
### Define subject identifiers ###
##################################
bblid=$1
dateid=$2
scanid=$3 

echo $bblid
echo $dateid
echo $scanid

subj="${bblid}"/"${dateid}"x"${scanid}"

## FreeSurfer subject directory
SUBJECTS_DIR=/data/joy/BBL/studies/pnc/processedData/structural/freesurfer53

## DSI Studio executable path
dsiBin="/share/apps/dsistudio/2016-01-25/bin/dsi_studio"

##############################################################
### Define directory structure and inputs for each subject ###
##############################################################
echo $bblid
echo $scanid
echo $dateid

roalfDir=/data/joy/BBL/studies/pnc/processedData/diffusion/pncDTI_2016_04/${bblid}/"${dateid}"x"${scanid}"/DTI_64
	
outDir=/data/joy/BBL/projects/pncBaumDti/GQI_test/${bblid}/"${dateid}"x"${scanid}"

logDir=${outDir}/logfiles

inputDir=${outDir}/input

dsi_outdir=${outDir}/dsiStudioRecon
	
baumDir=/data/joy/BBL/studies/pnc/processedData/diffusion/deterministic_20161201/"${bblid}"/"${dateid}"x"${scanid}"

mkdir -p ${dsi_outdir}
mkdir -p ${inputDir}/roi/Schaefer
mkdir -p ${inputDir}/roi/LausanneScale125

dilated_LausanneScale125_path="${baumDir}"/roi/LausanneScale125/${bblid}_"${dateid}"x"${scanid}"_LausanneScale125_dil4_dti.nii.gz

ln -s ${dilated_LausanneScale125_path} ${inputDir}/roi/LausanneScale125/

echo ${dilated_LausanneScale125_path}

#################################################################################
## Assign DTI image variables - Use eddy, motion, and distortion-corrected DWI ##
#################################################################################
# masked_dico_path=${roalfDir}/dico_corrected/"${bblid}"_"${dateid}"x"${scanid}"_dico_dico_mask.nii.gz	
# mask_path=${roalfDir}/dico_corrected/"${bblid}"_"${dateid}"x"${scanid}"_dico_dico_mask_mask.nii.gz

masked_dico_path=${baumDir}/input/"${bblid}"_"${dateid}"x"${scanid}"_masked_dico_dico.nii.gz
mask_path=${baumDir}/input/"${bblid}"_"${dateid}"x"${scanid}"_dtistd_2_mask.nii.gz

# Sym-link masked DTI file
ln -s ${masked_dico_path} ${inputDir}/
	
# Sym-link mask
ln -s ${mask_path} ${inputDir}/

################################################################
## Define subject-specific Rotated bvecs and other DTI inputs ##
################################################################
bvecs=${roalfDir}/raw_merged_dti/"${bblid}"_"${dateid}"x"${scanid}"_dti_merged_rotated.bvec
echo " "
echo "Subject-specific rotated bvecs file"
echo " "
echo ${bvecs}

ln -s ${bvecs} ${inputDir}/

## Bvals and acqparams 	are identical for all subjects ##
bvals=${roalfDir}/raw_merged_dti/"${bblid}"_"${dateid}"x"${scanid}"_dti_merged.bval 
echo " "
echo "bval file"
echo " "
echo ${bvals}
	
ln -s ${bvals} ${inputDir}/

indexfile=/data/joy/BBL/projects/pncReproc2015/diffusionResourceFiles/index_64.txt
acqparams=/data/joy/BBL/projects/pncReproc2015/diffusionResourceFiles/acqparams.txt 
	
##############################################
### Convert DWI input to DSI Studio format ###
##############################################
${dsiBin} --action=src --source="${masked_dico_path}" --bval="${bvals}" --bvec="${bvecs}" --output="${dsi_outdir}"/${bblid}_"${dateid}"x"${scanid}"_masked_dico_dico.src.gz

# Define DSI Studio source file for reconstruction and tractography
dti_source="${dsi_outdir}"/${bblid}_"${dateid}"x"${scanid}"_masked_dico_dico.src.gz

## DTI reconstruction
${dsiBin} --action=rec --source=${dti_source} --mask=${mask_path} --method=1 
	
## Define reconstruction file (fib)
dti_reconstruction=$(ls  ${dsi_outdir}/"${bblid}"_"${dateid}"x"${scanid}"_masked_dico_dico*dti*fib.gz )

## Export the qa0 and gfa maps from GQI fib file
${dsiBin} --action=exp --source=${dti_reconstruction} --export=fa0

########################
## GQI reconstruction ##
########################

################################################################################################
# Documentation: 
# http://dsi-studio.labsolver.org/Manual/Reconstruction#TOC-Generalized-Q-sampling-Imaging-GQI-
################################################################################################

${dsiBin} --action=rec --source=${dti_source} --mask=${mask_path} --method=4 --param0=1.25 --num_fiber=2 # --csf_calibration=1
	
## Define reconstruction file (fib)
gqi_reconstruction=$(ls ${dsi_outdir}/"${bblid}"_"${dateid}"x"${scanid}"_masked_dico_dico*gqi*fib.gz )
	
## Export the qa0 and gfa maps from GQI fib file
${dsiBin} --action=exp --source=${gqi_reconstruction} --export=fa0,gfa

####################################
## Run Deterministic Tractography ##
####################################

## Tractography output directory
tract_dir="${outDir}"/tractography
mkdir -p ${tract_dir}/connectivity
	
# DTI(FA-guided) tractography
${dsiBin} --action=trk --source="${dti_reconstruction}" --method=0 --fiber_count=1000000 --turning_angle=45 --min_length=10 --max_length=400 --output=${tract_dir}/${bblid}_"${dateid}"x"${scanid}"_dti_tractography.trk.gz --export="stat"

dti_tractography=${tract_dir}/${bblid}_"${dateid}"x"${scanid}"_dti_tractography.trk.gz

echo "DTI (FA-guided) Tractography output"
echo ""
echo ${dti_tractography}
echo ""

# GQI (QA-guided) tractography
${dsiBin} --action=trk --source="${gqi_reconstruction}" --method=0 --fiber_count=1000000 --turning_angle=45 --min_length=10 --max_length=400 --output=${tract_dir}/${bblid}_"${dateid}"x"${scanid}"_gqi_tractography.trk.gz --export="stat"

gqi_tractography=${tract_dir}/${bblid}_"${dateid}"x"${scanid}"_gqi_tractography.trk.gz

echo "GQI (QA-guided) Tractography output"
echo ""
echo ${gqi_tractography}
echo ""

###############################################
## Create Schaefer Parcellation in DTI space ##
###############################################

## Define inputs to ANTs call
PNC_template=/home/rciric/xcpAccelerator/xcpEngine/space/PNC/PNC-9375x9375x1.nii.gz
schaeferPNC_template=/data/joy/BBL/applications/xcpEngine/networks/SchaeferPNC.nii.gz
dti_refVol=/data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${bblid}/"${dateid}"x"${scanid}"/dti2xcp/${bblid}_"${dateid}"x"${scanid}"_referenceVolume.nii.gz
struct2seq_coreg=/data/joy/BBL/studies/pnc/processedData/diffusion/dti2xcp_201606230942/${bblid}/"${dateid}"x"${scanid}"/coreg/${bblid}_"${dateid}"x"${scanid}"_struct2seq.txt
antsDir=/data/joy/BBL/studies/pnc/processedData/structural/antsCorticalThickness/${bblid}/"${dateid}"x"${scanid}"
pncTemplate2subjAffine=${antsDir}/TemplateToSubject1GenericAffine.mat
pncTemplate2subjWarp=${antsDir}/TemplateToSubject0Warp.nii.gz


## Parcellations in dti space (output)
schaeferDti200_path=${inputDir}/roi/Schaefer/${bblid}_"${dateid}"x"${scanid}"_SchaeferPNC_200_dti.nii.gz
schaeferDti400_path=${inputDir}/roi/Schaefer/${bblid}_"${dateid}"x"${scanid}"_SchaeferPNC_400_dti.nii.gz
schaeferDti600_path=${inputDir}/roi/Schaefer/${bblid}_"${dateid}"x"${scanid}"_SchaeferPNC_600_dti.nii.gz
schaeferDti800_path=${inputDir}/roi/Schaefer/${bblid}_"${dateid}"x"${scanid}"_SchaeferPNC_800_dti.nii.gz
schaeferDti1000_path=${inputDir}/roi/Schaefer/${bblid}_"${dateid}"x"${scanid}"_SchaeferPNC_1000_dti.nii.gz

## Move MNI templates to PNC
schaeferDir=/data/joy/BBL/projects/pncBaumDti/Schaefer2018_LocalGlobal_Parcellation/CBIG-master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI
schaefer200_mni=${schaeferDir}/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_1mm.nii.gz
schaefer200_pnc=/data/joy/BBL/applications/xcpEngine/networks/SchaeferPNC_200.nii.gz
schaefer400_mni=${schaeferDir}/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz
schaefer400_pnc=/data/joy/BBL/applications/xcpEngine/networks/SchaeferPNC_400.nii.gz
schaefer600_mni=${schaeferDir}/Schaefer2018_600Parcels_7Networks_order_FSLMNI152_1mm.nii.gz
schaefer600_pnc=/data/joy/BBL/applications/xcpEngine/networks/SchaeferPNC_600.nii.gz
schaefer800_mni=${schaeferDir}/Schaefer2018_800Parcels_7Networks_order_FSLMNI152_1mm.nii.gz
schaefer800_pnc=/data/joy/BBL/applications/xcpEngine/networks/SchaeferPNC_800.nii.gz
schaefer1000_mni=${schaeferDir}/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_1mm.nii.gz
schaefer1000_pnc=/data/joy/BBL/applications/xcpEngine/networks/SchaeferPNC_1000.nii.gz

pncTransformDir=/home/rciric/xcpAccelerator/xcpEngine/space/PNC/PNC_transforms

## Move parcellations from MNI space to PNC template
antsApplyTransforms -e 3 -d 3 -i ${schaefer200_mni} -r ${PNC_template} -o ${schaefer200_pnc} -t ${pncTransformDir}/MNI-PNC_1Warp.nii.gz -t ${pncTransformDir}/MNI-PNC_0Affine.mat -n Multilabel

## Move Schaefer parcellation from PNC template to subject DTI space ##	
antsApplyTransforms -e 3 -d 3 -i ${schaefer200_pnc} -r ${dti_refVol} -o ${schaeferDti200_path} -t ${struct2seq_coreg} -t ${pncTemplate2subjAffine} -t ${pncTemplate2subjWarp} -n MultiLabel

## Create dilated Schaefer template in DTI space
dilated_schaeferDti200_path=${inputDir}/roi/Schaefer/${bblid}_"${dateid}"x"${scanid}"_SchaeferPNC_200_dti_dil2.nii.gz

## Dilate by 2 voxels in diffusion space (4mm)
ImageMath 3 ${dilated_schaeferDti200_path} GD ${schaeferDti200_path} 2 

##################################################
## Get volume of each ROI for each parcellation ##
##################################################
schaefer_vol_output=${inputDir}/roi/Schaefer/SchaeferPNC_200_dti_dil2_roiVol.txt

## Schaefer 200
for reg in {1..200}; do 
	echo ${reg}
	3dBrickStat -non-zero -count "${dilated_schaeferDti200_path}<${reg}>" 2>>/dev/null 1>> ${schaefer_vol_output}
done

lausanne_vol_output=${inputDir}/roi/LausanneScale125/LausanneScale125_dti_roiVol.txt

## Lausanne 234
for reg in {1..234}; do 
	echo ${reg}
	3dBrickStat -non-zero -count "${dilated_LausanneScale125_path}<${reg}>" 2>>/dev/null 1>> ${lausanne_vol_output}
done

###################################################################
## DSI Studio command to generate Schaefer connectivity matrices ##
###################################################################
"${dsiBin}" --action=ana --source="${gqi_reconstruction}" --tract="${gqi_tractography}" --connectivity="${dilated_schaeferDti200_path}" --connectivity_value=qa,count,mean_length --connectivity_type=end

## Rename and move matrices to proper output directory
Schaefer_qa_mat=$(ls ${dsi_outdir}/*SchaeferPNC*qa*.mat)
Schaefer_count_mat=$(ls ${dsi_outdir}/*SchaeferPNC*count*.mat)
Schaefer_length_mat=$(ls ${dsi_outdir}/*SchaeferPNC*length*.mat)

mv ${Schaefer_qa_mat} ${tract_dir}/connectivity/"${scanid}"_SchaeferPNC_200_qa_connectivity.mat
mv ${Schaefer_count_mat} ${tract_dir}/connectivity/"${scanid}"_SchaeferPNC_200_count_connectivity.mat
mv ${Schaefer_length_mat} ${tract_dir}/connectivity/"${scanid}"_SchaeferPNC_200_length_connectivity.mat

## Generate DTI-based (FA) matrices
"${dsiBin}" --action=ana --source="${dti_reconstruction}" --tract="${dti_tractography}" --connectivity="${dilated_schaeferDti200_path}" --connectivity_value=fa,count --connectivity_type=end

## Rename and move matrices to proper output directory
Schaefer_fa_mat=$(ls ${dsi_outdir}/*SchaeferPNC*fa*.mat)
Schaefer_count_mat=$(ls ${dsi_outdir}/*SchaeferPNC*count*.mat)

mv ${Schaefer_fa_mat} ${tract_dir}/connectivity/"${scanid}"_SchaeferPNC_200_dti_fa_connectivity.mat
mv ${Schaefer_count_mat} ${tract_dir}/connectivity/"${scanid}"_SchaeferPNC_200_dti_count_connectivity.mat

###################################################################
## DSI Studio command to generate Lausanne connectivity matrices ##
###################################################################
"${dsiBin}" --action=ana --source="${gqi_reconstruction}" --tract="${gqi_tractography}" --connectivity=${dilated_LausanneScale125_path} --connectivity_value=qa,count,mean_length --connectivity_type=end

LausanneScale125_qa_mat=$(ls ${dsi_outdir}/*LausanneScale125*qa*.mat)
LausanneScale125_count_mat=$(ls ${dsi_outdir}/*LausanneScale125*count*.mat)
LausanneScale125_length_mat=$(ls ${dsi_outdir}/*LausanneScale125*length*.mat)

mv ${LausanneScale125_qa_mat} ${tract_dir}/connectivity/"${scanid}"_LausanneScale125_qa_connectivity.mat
mv ${LausanneScale125_count_mat} ${tract_dir}/connectivity/"${scanid}"_LausanneScale125_count_connectivity.mat
mv ${LausanneScale125_length_mat} ${tract_dir}/connectivity/"${scanid}"_LausanneScale125_length_connectivity.mat

## Generate DTI-based (FA) matrices
"${dsiBin}" --action=ana --source="${dti_reconstruction}" --tract="${dti_tractography}" --connectivity=${dilated_LausanneScale125_path} --connectivity_value=fa,count --connectivity_type=end

LausanneScale125_fa_mat=$(ls ${dsi_outdir}/*LausanneScale125*fa*.mat)
LausanneScale125_count_mat=$(ls ${dsi_outdir}/*LausanneScale125*count*.mat)

mv ${LausanneScale125_fa_mat} ${tract_dir}/connectivity/"${scanid}"_LausanneScale125_dti_fa_connectivity.mat
mv ${LausanneScale125_count_mat} ${tract_dir}/connectivity/"${scanid}"_LausanneScale125_dti_count_connectivity.mat

# end
