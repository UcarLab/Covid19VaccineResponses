The Dogma-seq assay produced 4 groups of fastq files. For ATAC, GEX, HTOs and ADTs. These are shared through dbGaP (The database of Genotypes and Phenotypes)

ATAC and GEX are processed together using Cell Ranger ARC software for multiome.
You can refer to cellranger.sh file for processing individual libraries. 
Then you can refer to cellranger-aggr.sh file for aggreageting the individual outputs.

We used KITE to process ADTs (CITE-seq probes for cell surface proteins) and HTOs (Hash Tag Oligos, for multiplexing).
Please refer to the scripts inside KITE folder for details.

We experimented with different tools and settled on these, however feel free to try others.

