#!/bin/bash

#$ -N fastQC_demo
#$ -o /data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/fastqc/fastqc.out
#$ -e /data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/fastqc/fastqc.err
#$ -q pub8i
#$ -pe openmp 8-64
#$ -m beas
#$ -ckpt blcr

module load blcr
module load fastqc/0.11.7

# These files contain the full path and names of the data files to be analyzed.
DATA_FILENAMES_R1=/data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R1_filenames.txt
DATA_FILENAMES_R2=/data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R2_filenames.txt

QC_OUT_DIR=/data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/fastqc
QC_HTML_DIR=/data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/fastqc/fastqc_html
HTML="fastqc_html"

# Here we are performing a loop that will find each file in our data_filenames.txt file.
# Each file will be processed with the program "fastqc"
# "\" symbol indicates that more options for the program are on the next line 
# (--outdir) indicates the output directory for the result files
# (mv) moves the HTML result files to the new HTML result file directory

for FILE in $(cat ${DATA_FILENAMES_R1} ${DATA_FILENAMES_R2}); do
    fastqc $FILE \
    --outdir ${QC_OUT_DIR}
    
    mv ${QC_OUT_DIR}/*.html ${QC_HTML_DIR}
done

# Here we are compressing the HTML result file using the program tar
# -C flag prevents the parent directories from being included in the archive
# -csvf (c)reates archive, uses g(z)ip for compression, (v)erbosely shows the .tar file progress, (f)ilename appears next in the command
tar -C ${QC_OUT_DIR} -czvf ${HTML}.tar.gz ${HTML}


mv ${HTML}.tar.gz ${QC_HTML_DIR}/${HTML}.tar.gz
