# COVID-19 Vaccine Responses

This repository contains the code for our manuscript.


Single-cell data (annotations, RNA and ADT counts) is available on the [CZ CELLxGENE portal](https://cellxgene.cziscience.com/collections/5912a204-1095-4ea1-9ab1-362a638fe5a2).
The ADT (antibody-derived tag) protein counts are stored as an AnnData object in `adata.uns['ADT_protein_anndata']`.

Explore the data in our [interactive app](https://ucarlab.github.io/Covid19VaccineResponses).

Giray also shares some additional perspectives on his blog, [Good to Think About](https://goodtothinkabout.com). Opinions are his own.

## Repository contents

* `ATAC-seq/`: code for ATAC-seq analysis
* `FASTQ file processing/`: scripts and files for processing raw FASTQ files (see the README inside for details)
* `GEX-ADT/`: code for gene expression and surface protein data analysis
* `Reanalysis of PMID-34252919/`: Jupyter notebooks that reanalyze data from ["Systems vaccinology of the BNT162b2 mRNA vaccine in humans"](https://pubmed.ncbi.nlm.nih.gov/34252919/) (2021, PMID: 34252919)
