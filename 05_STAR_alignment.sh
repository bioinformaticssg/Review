#!/bin/bash

#$ -N PE_staralign_local
#$ -o /data/users/$USER/BioinformaticsSG/griffith_analysis_demo/alignments/staralign.out
#$ -e /data/users/$USER/BioinformaticsSG/griffith_analysis_demo/alignments/staralign.err
#$ -q pub8i
#$ -pe openmp 8
#$ -m beas
#$ -ckpt blcr

# Modified script orginally by A.Garibaldi (https://github.com/bioinformaticssg/Alignments/blob/master/staralign_transcriptcounts_local.sh)

set -euxo pipefail

module load blcr
module load STAR/2.5.2a
module load enthought_python/7.3.2
module load samtools/1.3

echo "You're aligning on $HOSTNAME"

# Number of threads. MUST match the qsub header openmp value you chose
P=8

# Path to all reference files
REF_DIR=/data/users/$USER/BioinformaticsSG/griffith_data/refs
# Path the reference genome
REF_FASTA=${REF_DIR}/chr22.ERCC92.fa
# Path to the annotation file
REF_ANNOTATION=${REF_DIR}/chr22.ERCC92.gtf
# Path to index files
INDEX_DIR=${REF_DIR}/index_files

# Directory for alignment result files
ALIGNMENTS_DIR=/data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/alignments

DATA_FILENAMES_R1=/data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R1_filenames.txt
DATA_FILENAMES_R2=/data/users/$USER/BioinformaticsSG/Review/griffith_analysis_demo/data/HBR_data_R2_filenames.txt
BASENAME="HBR"

sample_count=$(cat ${DATA_FILENAMES_R1} | wc -l)


#####################

# Iterate over each sample

for SAMPLE_N in $(seq ${sample_count}); do

    # Build the name of the files.
    R1=$(head -n ${SAMPLE_N} ${DATA_FILENAMES_R1} | tail -n 1)
    R2=$(head -n ${SAMPLE_N} ${DATA_FILENAMES_R2} | tail -n 1)

    # Run the aligner.
	echo "Aligning ${R1} and ${R2}"

	STAR \
	--chimSegmentMin 15 \
	--chimJunctionOverhangMin 15 \
	--outFilterMismatchNmax 3 \
	--alignEndsType Local \
	--runThreadN $P \
	--outFilterMultimapNmax 1 \
	--outBAMcompression 10 \
	--outBAMsortingThreadN $P \
	--quantMode GeneCounts  TranscriptomeSAM  \
	--quantTranscriptomeBAMcompression 10 \
	--outSAMtype BAM SortedByCoordinate \
	--alignSJDBoverhangMin 6 \
	--alignIntronMax 300000 \
	--genomeDir $INDEX_DIR \
	--sjdbGTFfile $REF_ANNOTATION \
	--outFileNamePrefix ${ALIGNMENTS_DIR}/${BASENAME}_${SAMPLE_N}.star \
	--readFilesCommand zcat \
	--readFilesIn ${R1} ${R2}

	echo "Alignment of ${R1} and ${R2} complete"
done

exit
########################################

#STAR --chimSegmentMin 15 --chimJunctionOverhangMin 15 \ #allows for circRNA detection
#--outFilterMismatchNmax 3 \ #maximum mismatches allowed
#--alignEndsType Local  --runThreadN $P \ #type of alignment and number of threads.Local is standard local alignment with soft-clipping allowed
#--quantMode GeneCounts  TranscriptomeSAM  --quantTranscriptomeBAMcompression 10 \ #performs gene counts like HTSEQ union mode. Just exons counts. max compression
#--outFilterMultimapNmax 1 \ # default is 10.setting this to 1 will limit the bam file to contain only TRULY uniquely mapped reads. Stats will still contain multimap info, but the definition of a multipmapper more correct by setting this to 1 instead of 10 $
#--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
#--outBAMcompression 10 --outBAMsortingThreadN $P \ #these features increase compression and compression speed by multithreading more than the default 6 threads
#--alignSJDBoverhangMin 6 --alignIntronMax 300000 \ #default intron is 0. if this is 0 it turns off splicing detection.  default overhang is 5

#NOTES:
#--outSAMstrandField intronMotif needed if non stranded data
#.run this if error produced. input error amount --limitBAMsortRAM 30943606211
#note if you run out of ram, tactic is to take out BAM sorted option to yield unsorted then do it separately with samtools