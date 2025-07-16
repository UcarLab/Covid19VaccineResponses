#!/bin/sh

#SBATCH --job-name=agg
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=70
#SBATCH --mem=600GB
#SBATCH -q long
#SBATCH -p compute

cd /projects/ucar-lab/giray/DOGMA-seq/Dogma-2

REF_FILE='refdata-cellranger-arc-GRCh38-2020-A-2.0.0'
OUTPUT_FOLDER="aggregated_samples"
AGGR_FILE="aggregation.csv"

echo "The reference file is: $REF_FILE"
echo "The aggregation file (the csv) is: $AGGR_FILE"
echo "The output folder is: $OUTPUT_FOLDER"

cellranger-arc-2.0.1/cellranger-arc aggr \
                                            --id="$OUTPUT_FOLDER" \
                                            --reference="$REF_FILE" \
                                            --csv="$AGGR_FILE" \
                                            --disable-ui \
                                            --normalize=none \
                                            --nosecondary

