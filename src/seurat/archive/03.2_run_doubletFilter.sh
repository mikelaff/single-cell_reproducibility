#!/bin/bash

# files dir
filesDir="/proj/steinlab/projects/IVIV_scRNA/single-cell_reproducibility/results/seurat/PGP1_QC_split_by_sample_seurat_objects/"

###################################################################
## run splitp to change the hexamer barcode into oligodT barcode ##
###################################################################

R2_splitp="${dirName}_R2_barcodeCollapsed.fastq.gz"

if [ ! -f "${R2_splitp}" ]; then

	splitp \
	-r ${R2} \
	-b ${barcodeFile} \
	-s 79 \
	-e 86 \
	-o | gzip > ${R2_splitp}
fi

###################################################################
##                 fastp : removing low-quality reads            ##
###################################################################

module add fastp/0.23.2

R1_fastp="${dirName}_R1_fastp.fastq.gz"
R2_fastp="${dirName}_R2_fastp.fastq.gz"

fastp \
--disable_adapter_trimming \
-i ${R1} \
-I ${R2_splitp} \
-o ${R1_fastp} \
-O ${R2_fastp} \
-w ${numThread} \
-j ${dirName} \
-h ${dirName} \
-q 20 \
-n 1  

###################################################################
##################        Run alevin-fry       ####################
###################################################################

module add salmon/2.0.1
module add alevin-fry/0.5.1

	echo "running salmon alevin with 2[1-100]"

salmon alevin \
        -i ${indexDir} \
        -l A \
        -1 ${R1_fastp} \
        -2 ${R2_fastp} \
	--splitseqV2 \
	-p ${numThread} \
	-o ${dirName}_alevin \
	--sketch

alevin-fry generate-permit-list \
-d both \
-i ${dirName}_alevin \
--output-dir ${dirName}_alevin/permit_list \
--unfiltered-pl ${barcodeComb} \
--min-reads 100

alevin-fry collate \
-r ${dirName}_alevin \
-i ${dirName}_alevin/permit_list \
-t ${numThread} \
--compress

alevin-fry quant \
-m ${refTranscriptome} \
-i ${dirName}_alevin/permit_list \
-t ${numThread} \
-o ${dirName}_alevin/countMatrix \
-r cr-like \
--use-mtx





