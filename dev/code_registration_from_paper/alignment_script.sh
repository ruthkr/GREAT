#!/bin/bash

#SBATCH -p jic-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of tasks
#SBATCH -t 1-0:00 #1-0:00 # time (D-H:m)
#SBATCH -o ./slurm_output/slurm.%j.out # STDOUT
#SBATCH -e ./slurm_output/slurm.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications by email in the event of these events
#SBATCH --mail-user=alexander.calderwood@jic.ac.uk # emails sent to address
#SBATCH --mem 64000 #64288 # memory pool for ALL cores
#SBATCH --localscratch=ssd:120 # request SSD harddrive space (Gb)

# Tool versions used:
# source fastqc-0.11.3
# source hisat-2.0.4
# source samtools-1.5
# source stringtie-1.2.2
# source bbtools-37.68

set -e #should report as failiure + quit if any step fails
set -u
set -x

SAMPLE_ID=$1 #E.G. A1
DS_ID=$2 # E.G. ds_17_11_13
QC_DIR=$3 # E.G. ../QC
RAW_DIR=$4 # E.G. ../raw_data
#INTERMEDIATE_DIR=$5 # E.G. ../intermediate_data
INTERMEDIATE_DIR=../intermediate_data #${SLURM_LOCAL_SCRATCH} # use local SSD on machine!
FINAL_DIR=$6 # E.G. ../final_data
REF_URI=$7 # E.G. ../reference_data/AC_genome
REF_GFF=$8

echo "sample: $SAMPLE_ID"

# rm -rf ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}*
mkdir -p ${INTERMEDIATE_DIR}/${DS_ID}

#init ial raw read QC
echo "doing raw QC..."
srun fastqc -o ${QC_DIR}/raw_reads/${DS_ID}/FastQC_reports/ \
${RAW_DIR}/${DS_ID}/data/${SAMPLE_ID}_1.fq.gz \
${RAW_DIR}/${DS_ID}/data/${SAMPLE_ID}_2.fq.gz

#trim first 15bp from all reads.
#Produces P files (reads which still paired after trimming.)
# and U files, reads where only 1 of pair remains after trimming.
# all U files are empty, as all reads remain paired after trimming 15bp, therefore
# can be discared.
# v0.33
echo "trimming first 15 bp..."
srun java -jar "/nbi/software/testing/trimmomatic/0.33/x86_64/bin/trimmomatic-0.33.jar" PE \
-phred64 \
-basein ${RAW_DIR}/${DS_ID}/data/${SAMPLE_ID}_1.fq.gz \
-baseout ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_trimmed.fq.gz \
HEADCROP:15

# clear up empty _U files created.
rm ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_trimmed_1U.fq.gz
rm ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_trimmed_2U.fq.gz

echo "bbduk adaptor trimming..."
# Adaptor trimming. Adaptors will only ever be on 3' end. ktrim=r is normal mode
# for adaptor trimming.
# k sepcified kmer size to use ( at most, length of adaptors)
# mink allows use of shorter kmers at the ends of the read.
# hdist is hamming dist mismatch allowed for trimming. Hdist2 is like hdist, but
# is only applied to kmeers shorter than k (down to min k)
# tbo -> also trim adaptors based on pair overlap detection and tpe - trim both
# #reads of paired to same length - recommneded by the manual
# v
srun bbduk.sh \
in1=${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_trimmed_1P.fq.gz \
in2=${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_trimmed_2P.fq.gz \
out1=${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_more_trimmed_1P.fq.gz \
out2=${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_more_trimmed_2P.fq.gz \
ref=Novogene_Illumina_adaptors.fa,truseq_rna.fa.gz,truseq.fa.gz,adapters.fa \
ktrim=r \
k=23 \
mink=11 \
hdist=1 \
hdist2=0 \
tbo \
tpe

# QC on adaptor trimmed reads
echo "doing trimmed QC..."
srun fastqc -o ${QC_DIR}/raw_reads/${DS_ID}/FastQC_reports_trimmed/ \
${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_more_trimmed_1P.fq.gz \
${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_more_trimmed_2P.fq.gz

#hisat align
echo "doing hisat2 align with strandedness..."
# nb as libraries produced using "NEB next ultra directional library kit",
# have stranded reads.
# --dta : output includes XS tag to indicate genomic strand producing.
#-x sould point to the hisat index base name
srun hisat2 --dta \
-x ${REF_URI} \
-1 ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_more_trimmed_1P.fq.gz \
-2 ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}_more_trimmed_2P.fq.gz \
-S ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sam \
--rna-strandness RF

# sort and convert to bam
echo "sorting, and converting to bam..."
srun samtools sort -o ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.bam \
-O bam \
${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sam
#index
echo "indexing..."
srun samtools index ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.bam
# generate flagstat alignment stats, and to QC
srun samtools flagstat ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.bam > \
${QC_DIR}/read_alignment/${DS_ID}/raw_aligned/${SAMPLE_ID}.txt
#cleanup
srun rm ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sam

#filter to keep only properly pairing
echo "filtering aligned reads..."
# see https://broadinstitute.github.io/picard/explain-flags.html for decoding
# samtools filtering flags.
#-F gives the things to skip
#-f gives the things to keep
#-f 2 only include paired reads, with proper orientation.
#filter
srun samtools view -bh -f 2 \
-o ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.filtered.bam \
${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.bam
#sort
srun samtools sort \
-o ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.filtered.bam \
${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.filtered.bam
#index
srun samtools index \
${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.filtered.bam
# generate flagstat alignment stats, and to QC
srun samtools flagstat \
${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.filtered.bam > \
${QC_DIR}/read_alignment/${DS_ID}/filtered/${SAMPLE_ID}.txt
#cleanup
#srun rm ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.bam
#srun rm ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.bam.bai
#srun rm ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.filtered.bam


#assemble and quantify
echo "assembling with stringtie..."
#assemble & quantify reads
# -o <out.gtf> where write assembled reads to
# -p number of threads
# -G use the reference annotation file (.gtf or .gff3 format). If provided,
# stringtie will prefer to use these "known" genes. It will also produce
# additional transcrtipts to account for data that aren't covered by the
#annotation.
# -l <label> sets label as the the prefix for the output transcripts (default STRG)
# -A <gene_abund.tab> gene abundances reprted to where specified
# -B enables output of Ballgown input table files (*.ctab). Requires -G option.
# -e limits processing of read alignments to only transcripts given in -G file.
srun stringtie ${INTERMEDIATE_DIR}/${DS_ID}/${SAMPLE_ID}.sorted.filtered.bam \
-e \
-B \
-A ${FINAL_DIR}/${DS_ID}/ballgown/ERR_${DS_ID}_${SAMPLE_ID}/${DS_ID}_${SAMPLE_ID}_gene_abund.tab \
-G ${REF_GFF} \
-o ${FINAL_DIR}/${DS_ID}/ballgown/ERR_${DS_ID}_${SAMPLE_ID}/${DS_ID}_${SAMPLE_ID}_transcriptome.gtf \
-l ${SAMPLE_ID}
