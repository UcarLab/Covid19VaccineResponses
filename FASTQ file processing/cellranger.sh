#!/bin/sh

#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=250GB
#SBATCH --qos=batch

REF_FILE='refdata-cellranger-arc-GRCh38-2020-A-2.0.0'

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
    '2022-03-16-1-1'
    '2022-03-16-1-2'
    '2022-03-16-1-3'
    '2022-03-16-1-4'
    '2022-03-16-1-5'
    '2022-03-16-1-6'
    '2022-03-16-1-7'
    '2022-03-16-1-8'
    '2022-03-16-2-1'
    '2022-03-16-2-2'
    '2022-03-16-2-3'
)



for LIBRARIES_ID in ${ids[@]}; do

    OUTPUT_FOLDER="$LIBRARIES_ID-out"
    LIBRARIES_FILE="library_csv_files/$LIBRARIES_ID-libraries.csv"
    
    echo "Starting $LIBRARIES_ID"
    echo "The reference file is: $REF_FILE"
    echo "The libraries are noted in: $LIBRARIES_FILE"
    echo "The output folder is $OUTPUT_FOLDER"
    
    cellranger-arc-2.0.1/cellranger-arc count \
                                                --id="$OUTPUT_FOLDER" \
                                                --reference="$REF_FILE" \
                                                --libraries="$LIBRARIES_FILE" \
                                                --disable-ui
    
    
    echo 'Done.'

done
