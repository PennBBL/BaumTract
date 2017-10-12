#!/bin/sh

#################################################################
### Define subject list in 'bblid_scanid_dateidofscan' format ###
#################################################################
subject_list=$(cat /data/joy/BBL/projects/pncBaumDti/subject_lists/LTN_pncReproc_baumDTI_subjList_n879.txt)

# /data/jag/gbaum/PNC/10_subj_test/DTI_R21_TXFR/full_80010_subjID.txt)

# /data/joy/BBL/projects/pncBaumDti/Motion_paper/subjects/subject_lists/n949_dtiQA_b0_FSexclude_LTNv2_subjList.txt)
for name in ${subject_list}; do

	bblid=$(basename ${name} | cut -d_ -f1)
	scanid=$(basename ${name} | cut -d_ -f2)
	dateid=$(basename ${name} | cut -d_ -f3)

	echo $bblid
	echo $scanid
	echo $dateid

	subj="${bblid}"/"${dateid}"x"${scanid}"

	#####################################
	### Define input and output paths ###
	#####################################
	roalfDir=/data/joy/BBL/studies/pnc/processedData/diffusion/pncDTI_2016_04/${bblid}/"${dateid}"x"${scanid}"/DTI_64
	
	outDir=/data/joy/BBL/projects/pncBaumDti/GQI_test/${bblid}/"${dateid}"x"${scanid}"

	logDir=${outDir}/logfiles

	inputDir=${outDir}/input

	dsi_outdir=${outDir}/dsiStudioRecon
		
	baumDir=/data/joy/BBL/studies/pnc/processedData/diffusion/deterministic_20161201/"${bblid}"/"${dateid}"x"${scanid}"

	mkdir -p ${outDir}/logfiles

 	var0="pushd /data/joy/BBL/projects/pncBaumDti/GQI_test/scripts; ./process_GQI.sh ${bblid} ${dateid} ${scanid}; popd"

	# Remove old scripts/logs
	rm ${logDir}/run_processGQI_detPipe_"${bblid}"_"${dateid}"x"${scanid}".*

	# Write command to subject-specific executable script
	# echo -e "${var0}"
	echo -e "${var0}" >> ${logDir}/run_processGQI_detPipe_"${bblid}"_"${dateid}"x"${scanid}".sh

	subject_script=${logDir}/run_processGQI_detPipe_"${bblid}"_"${dateid}"x"${scanid}".sh
	
	chmod 775 ${subject_script}
 	
	## Execute qsub job for probtrackx2 runs for each subject 
	qsub -q all.q,basic.q -wd ${logDir} -l h_vmem=6G,s_vmem=5G ${subject_script}

done
