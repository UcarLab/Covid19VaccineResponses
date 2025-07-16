#!/bin/sh

#SBATCH --job-name=adt
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250GB
#SBATCH --qos=batch


ids=(
    '2023-02-20-1-1'
    '2023-02-20-1-2'
    '2023-02-20-1-3'
    '2023-02-20-1-5'
    '2023-02-20-1-6'
    '2023-02-20-1-7'
    '2023-02-21-1-4'
    '2023-02-21-1-5'
    '2023-02-21-1-6'
    '2023-02-21-1-7'
    '2023-03-01-1-3'
    '2023-03-01-2-2'
    '2023-03-02-1-1'
    '2023-03-02-1-2'
    '2023-03-02-1-3'
    '2023-02-21-1-2'
    '2023-03-01-1-4'
    '2023-02-21-1-3'
)

fastqs_folder="RAWDATA"
kallisto index -i FeaturesMismatch.idx -k 15 ./FeaturesMismatch.fa

white_list="/projects/ucar-lab/giray/DOGMA-seq/737K-arc-v1.txt"


for LIBRARIES_ID in ${ids[@]}; do

    mkdir $LIBRARIES_ID

    kallisto bus -i FeaturesMismatch.idx -o $LIBRARIES_ID -x 10xv3 -t $SLURM_CPUS_PER_TASK "$fastqs_folder"/FeatureHTO-"$LIBRARIES_ID"*.fastq.gz
    
    bustools correct -w "$white_list"  "$LIBRARIES_ID"/output.bus -o "$LIBRARIES_ID"/output_corrected.bus
    
    bustools sort -t $SLURM_CPUS_PER_TASK -o "$LIBRARIES_ID"/output_sorted.bus "$LIBRARIES_ID"/output_corrected.bus
    
    bustools count -o "$LIBRARIES_ID"/featurecounts --genecounts -g ./FeaturesMismatch.t2g -e "$LIBRARIES_ID"/matrix.ec -t "$LIBRARIES_ID"/transcripts.txt "$LIBRARIES_ID"/output_sorted.bus
    

done