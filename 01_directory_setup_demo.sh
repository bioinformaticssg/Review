#!/bin/bash

# set tells the program to stop if it runs into any issues, the flag descriptions are as follows:
# -e Exit immediately when a command fails.
# -u Treats unset/unbound variables as an error and exits immediately.
# -x Prints each command before executing it -- this is helpful for debugging.
# -o Sets the exit code (0 is successful) to the that of the rightmost command...continued on next line
# -o If the first command in a pipeline failed it will be carried through to the end of the pipeline and still exit instead of continuing through the script

# Change experiment name
export EXP_NAME=griffith_analysis_demo
# Change root directory
export ROOT_DIR=/data/users/$USER/BioinformaticsSG
# Change sub directory
export SUB_DIR=${ROOT_DIR}/Review
# Change location of reference genome
export REF_DIR=${ROOT_DIR}/griffith_data/refs
# Change to indexed reference genome
export GENOME="22"

############################################################################################

# These variables identify the paths for the genome reference and gene features files 
export REF=${REF_DIR}/${GENOME}
export REF_FA=${REF}.fa
export REF_GTF=${REF}.gtf

# These variables identify the paths for the directories for your project result files
# This is your project directory
export EXP_DIR=${SUB_DIR}/${EXP_NAME}
# This is your scripts directory
export SCRIPT_DIR=${EXP_DIR}/scripts
# This is the directory for your data
export DATA_DIR=${EXP_DIR}/data
# This is the directory for your alignment result files
export ALIGNMENTS_DIR=${EXP_DIR}/alignments
# This is the directory for your abundance estimation result files
export COUNTS_DIR=${EXP_DIR}/counts
# This is the directory for your differential expression (DE) analysis result files
export DE_DIR=${EXP_DIR}/DE_analysis

# This is the directory for your FastQC result files
export QC_OUT_DIR=${EXP_DIR}/fastqc_out
# This is the subdirectory for your HTML FastQC result files
export QC_HTML_DIR=${QC_OUT_DIR}/fastqc_html

# This is the direcotry for your trimmed data anlaysis
export TRIM_DIR=${EXP_DIR}/trimmed_analysis
# This is the directory for your trimmed data
export TRIM_DATA_DIR=${TRIM_DIR}/trimmed_data
# This is the directory for your FastQC result files for your trimmed data
export TRIM_QC_DIR=${TRIM_DIR}/trimmed_fastqc
# This is the subdirectory for your HTML FastQC result files for your trimmed data
export TRIM_QC_HTML_DIR=${TRIM_QC_DIR}/trim_fastqc_html
# This is the directory for your alignment result files from your trimmed data
export TRIM_ALIGNMENTS_DIR=${EXP_DIR}/trimmed_alignments
# This is the directory for your abundance estimated result files from your trimmed data
export TRIM_COUNTS_DIR=${EXP_DIR}/trimmed_counts

# Below is an example of setting a path to file that you may need for the program you are using.
# These files are used to remove the adapters from the sequences. 
# export SE_ADAPTER_DIR=/home/stacey/miniconda3/pkgs/trimmomatic-0.36-5/share/trimmomatic-0.36-5/adapters/TruSeq2-SE.fa
# export PE_ADAPTER_DIR=/home/stacey/miniconda3/pkgs/trimmomatic-0.36-5/share/trimmomatic-0.36-5/adapters/TruSeq2-PE.fa

# Here we make the directories that may not already exist.
# -p is a flag that tells the the mkdir program to ignore any warnings such as the directory already existing.
mkdir -p ${SCRIPT_DIR}
mkdir -p ${DATA_DIR}
mkdir -p ${ALIGNMENTS_DIR}
mkdir -p ${COUNTS_DIR}
mkdir -p ${DE_DIR}

mkdir -p ${QC_OUT_DIR}
mkdir -p ${QC_HTML_DIR}

mkdir -p ${TRIM_DIR}
mkdir -p ${TRIM_DATA_DIR}
mkdir -p ${TRIM_QC_DIR}
mkdir -p ${TRIM_QC_HTML_DIR}
mkdir -p ${TRIM_ALIGNMENTS_DIR}
mkdir -p ${TRIM_COUNTS_DIR}