This repository contains instructions and code to run the Variant-to-Phenotype (V2P) variant effect predictor. For many users, accessing V2P 
via our website, v2p.com, will be the best option. There, users may process up to 100,000 variants at a time. For advanced users or those who need 
higher throughput, offline installation may be a better option.

# Installing V2P

## Download the V2P repository

```git clone ```
Set the V2P working directory to the location that you cloned the repository.
```export V2P_DIR=/path/to/dir```
You can add this to your bash config if you would like it to persist.
```echo "export V2P_DIR=/path/to/dir" >> ~/.bashrc```

## Install dependencies

V2P's dependencies may be installed with conda. 
```conda env create --file v2p.yaml```
Activate the conda environment.
```conda activate v2p```
An additional python dependency must be installed after creating the conda environment.
this only needs to be done once. 
```
cd bin/scikit-multilearn
pip install .
```
V2P also requires an installation of docker to be available on the system.

### Download V2P data

V2P relies on a large set of features to derive its predictions. These features are collected from Ensembl's VEP and other sources.
In total, downloading the feature data requires ~590GB of free disk space. 
```
wget -O ${V2P_DIR}/.vep "link to vep data"
wget -O ${V2P_DIR}/hpo.db "link to hpo db"
# This can be place anywhere on your system. You must provide the path to this data
# when running V2P
wget "link to cadd dbs"
```

## (Optional) Download precomputed variants

Precomputed predictions for all posible single nucleotide variants and gnomAD indels are available for download. 
Installation will result in greatly increased speed for most variant sets. However, there is an additional 
disk space requirement (~806GB).
```
wget "link to precomputed"
```

## Running V2P

The input to V2P is a VCF file containing the variants you wish to score in hg38 coordinates. 
The VCF file must contain "#CHROM   POS ID  REF ALT QUAL    FILTER  INFO" in the header. 
Other header lines are not required.

To run V2P provide the path to the VCF, a path where the output will be stored, and the path to the downloaded features.
You may optionally pass a number of CPUs to be used (default=1) and, if scoring variants from a single gene, you may provide 
a gene name to be prioritized by VEP. Run ```bash get_predictions.sh -h``` for a description of the available parameters.

Example
```
bash get_predictions.sh -i /path/to/input -o /path/to/output -a /path/to/annotations -c 5
```