#!/bin/bash

#$ -N trim_qc_demo          # name of the job
#$ -o /data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/trim_data/trim_qc_demo.out   # contains what would normally be printed to stdout
#$ -e /data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/trim_data/trim_qc_demo.err   # file name to print standard error messages to. 
#$ -q pub8i                 # request cores from the free64, som, asom queues.
#$ -pe openmp 8             # request parallel environment. You can include a minimum and maximum core count.
#$ -m beas                  # send you email of job status (b)egin, (e)rror, (a)bort, (s)uspend
#$ -ckpt blcr               # (c)heckpoint: writes a snapshot of a process to disk, (r)estarts the process after the checkpoint is c$

set -euxo pipefail

module load blcr
module load fastqc/0.11.7

# These files contain the full path and names of the data files to be analyzed.
DATA_FILENAMES_R1=/data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R1_filenames.txt
DATA_FILENAMES_R2=/data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R2_filenames.txt
BASENAME="HBR"

DIR=/data/users/$USER/BioinformaticsSG/Reviewgriffith_analysis_demo

TRIM_DATA_DIR=${DIR}/trim_data
TRIM_QC_DIR=${TRIM_DATA_DIR}/trim_fastqc_out
TRIM_QC_HTML_DIR=${TRIM_QC_DIR}/trim_fastqc_html
HTML="trim_fastqc_html"

TRIMMOMATIC_DIR=/data/apps/trimmomatic/0.35/trimmomatic-0.35.jar 


# TRIMMOMATIC for paired end samples

for SAMPLE_N in $(seq 1 3); do

        # Build the name of the files.
        R1=$(head -n ${SAMPLE_N} ${DATA_FILENAMES_R1} | tail -n 1)
        R2=$(head -n ${SAMPLE_N} ${DATA_FILENAMES_R2} | tail -n 1)
        OUTPUT=${TRIM_DATA_DIR}/${BASENAME}_${SAMPLE_N}.fq.gz

        TRIMMER="HEADCROP:1"

        java -jar ${TRIMMOMATIC_DIR} \
        PE \
        -threads 8 \
        ${R1} ${R2} \
        -baseout ${OUTPUT} \
        ${TRIMMER}

     done
done

# FastQC on  trimmed data - paired data files only as indicated by the '\*P.\*'

for SAMPLE in `find ${TRIM_DATA_DIR} -name \*P.\*`; do
    fastqc ${SAMPLE} \
    --outdir ${TRIM_QC_DIR}

    # I am moving all the html files to a new directory
    mv ${TRIM_QC_DIR}/*.html ${TRIM_QC_HTML_DIR}
done

# I am compressing the directory containing all the html files into one file
tar -C ${TRIM_QC_DIR} -czvf ${HTML}.tar.gz ${TRIM_QC_HTML_DIR} 

mv ${HTML}.tar.gz ${TRIM_QC_HTML_DIR}

# Some notes on the trimmer setting:

# quality: Specifies the minimum quality required to keep a base.

# LEADING:<quality>
# Removes leading low quality or N bases (below quality 3)

# TRAILING:<quality>
# Remove trailing low quality or N bases (below quality 1)

# SLIDINGWINDOW:<windowSize>:<requiredQuality>
# Scan the read with a 4-base wide sliding window (windowSize), cutting when the average quality per base drops below 15 (requiredQuality)

# MINLEN:<length>
# Drop reads which are less than 36 bases long after these steps