#!/bin/bash

#$ -N star_build
#$ -o /data/users/$USER/BioinformaticsSG/griffith_data/refs/index_files/star_build.out 
#$ -e /data/users/$USER/BioinformaticsSG/griffith_data/refs/index_files/star_build.err 
#$ -q pub8i
#$ -pe openmp 8
#$ -m beas
#$ -ckpt blcr

# Modified script orginally by A.Garibaldi(/data/users/$USER/BioinformaticsSG/Alignments/star_build.sh)

set -euxo pipefail

module load blcr
module load STAR/2.5.2a
module load enthought_python/7.3.2
module load samtools/1.3

REF_DIR=/data/users/$USER/BioinformaticsSG/griffith_data/refs

# Path the reference genome
REF_FASTA=${REF_DIR}/chr22.ERCC92.fa
# Path to the annotation file
REF_ANNOTATION=${REF_DIR}/chr22.ERCC92.gtf
# Path to index files
INDEX_DIR=${REF_DIR}/index_files

O=99 #this overhang ideally should be ReadLength-1.
P=8 #threads

mkdir -p ${INDEX_DIR}

# Use the alignment program STAR to build the genome indices
STAR \
--runMode genomeGenerate \
 --genomeDir $OUTPUT_DIR \
 --genomeFastaFiles $REF_FASTA \
 --runThreadN $P \
 --sjdbGTFfile $REF_ANNOTATION \
 --sjdbOverhang $O