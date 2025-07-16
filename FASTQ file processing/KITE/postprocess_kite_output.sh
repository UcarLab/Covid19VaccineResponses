#!/bin/bash

set -e
## This is needed for the seurat function Read10X to be able to read the data.

## Side note:
## scanpy does not like the tranpose part and does not have a functioning read 10x function for this.
## One has to read the 2 txts and the mtx file separately using scanpy pandas etc.

[[ -z "$1" ]] && { echo "The (input/output) folder parameter is needed." ; exit 1; }

folder="$1"

echo "$folder"


python transpose_rename_kite_output.py "$folder/featurecounts.mtx" "$folder/matrix.mtx"
gzip -c "$folder/matrix.mtx" > "$folder/matrix.mtx.gz"
gzip -c "$folder/featurecounts.genes.txt" > "$folder/features.tsv.gz"
gzip -c "$folder/featurecounts.barcodes.txt" > "$folder/barcodes.tsv.gz"
