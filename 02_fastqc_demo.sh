#!/bin/bash

#$ -N fastQC_demo           # name of the job
#$ -o /data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/fastqc_out/fastqc.out  # contains what would normally be printed to stdout (the$
#$ -e /data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/fastqc_out/fastqc.err  # file name to print standard error messages to. These m$
#$ -q free64,som,asom       # request cores from the free64, som, asom queues.
#$ -pe openmp 8-64          # request parallel environment. You can include a minimum and maximum core count.
#$ -m beas                  # send you email of job status (b)egin, (e)rror, (a)bort, (s)uspend
#$ -ckpt blcr               # (c)heckpoint: writes a snapshot of a process to disk, (r)estarts the process after the checkpoint is c$

module load blcr
module load fastqc/0.11.7

DATA_FILENAMES_R1=/data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R1_filenames.txt
DATA_FILENAMES_R2=/data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R2_filenames.txt
QC_OUT_DIR=/data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/fastqc_out
QC_HTML_DIR=/data/users/sborrego/BioinformaticsSG/Review/griffith_analysis_demo/fastqc_out/fastqc_html
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
