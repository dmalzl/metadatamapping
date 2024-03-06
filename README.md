# metadatamapping
![pypi](https://img.shields.io/badge/pypi-v1.2.0-blue)
![python-version](https://img.shields.io/badge/Python->=3.10-blue)
![stable-version](https://img.shields.io/badge/version-1.2.0-blue)

A python library to fetch metadata from NCBI and MetaSRA for a list of NCBI accessions and data extraction from ARCHS4

## Installation
Simply run the following
```commandline
pip install metadatamapping
```
or clone the repository 
```commandline
git clone git@github.com:dmalzl/metadatamapping.git
```
and run
```commandline
cd metadatamapping
pip install .
```
you should then be able to import the package as usual

## Example usage
NCBI sample, experiment, biosample or geo accessions can be mapped to SRA uids using the `map_accessions_to_srauids` function from the `metadata` module of the package. The call to the function as shown below invokes two processes that concurrently fetch the SRA UIDs for the accessions in batches and write the results to the outputfile `"/path/to/outputfile"`.
```python
from metadatamapping import metadata
sra_uids = metadata.map_accessions_to_srauids(
    accessions,
    "/path/to/outputfile",
    n_processes = 2
)
```
The resulting SRA UIDs can then either be used to retrieve all associated accessions from the SRA with the `srauids_to_accessions` function from the `metadata` module like so
```python
from metadatamapping import metadata
ncbi_accessions = metadata.srauids_to_accessions(
    sra_uids
)
```
or link them to BioSample UIDs and then retrieve the associated metadata with the `link_sra_to_biosample` function from the `link` module and the `biosampleuids_to_metadata` function from the metadata module
```python
from metadatamapping import metadata, link
srauids_to_biosampleuids = link.link_sra_to_biosample(
    sra_uids.uid
)
biosample_metadata = metadata.biosample_uids_to_metadata(
    srauids_to_biosampleuids.biosample
)
```
Finally we can retrieve normalized metadata for the samples from [MetaSRA](https://metasra.biostat.wisc.edu/) using the `metasra_from_study_id` function of the `metadata` module (note that this database might not contain data for all your samples so the function may only returns normalized metadata for some of your samples)
```python
from metadatamapping import metadata
metasra_metadata = metadata.metasra_from_study_id(
    ncbi_accessions.study.unique()
)
```
Additionally, the package provides an interface for parsing the [ARCHS4 HDF5 format](https://maayanlab.cloud/archs4/help.html) which is located in the `archs4` module and handles parsing of associated metadata with the `get_filtered_sample_metadata` function as well as extraction of expression data in the [`AnnData`](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html) format with the `samples` function
```python
archs4_file = "/path/to/archs4.h5"
retain_keys = [
    'geo_accession', 'characteristics_ch1', 'molecule_ch1', 'readsaligned', 'relation', 
    'series_id', 'singlecellprobability', 'source_name_ch1', 'title'
]
archs4_metadata = archs4.get_filtered_sample_metadata(
    archs4_file,
    retain_keys
)
archs4_adata = archs4.samples(
    archs4_file,
    dataframe_indexed_by_geo_accessions,
    n_processes = 2
)
```
For a full demonstration of usage please refer to the [`Snakefile`](https://github.com/dmalzl/metadatamapping/blob/main/examples/Snakefile) in the [`examples`](https://github.com/dmalzl/metadatamapping/tree/main/examples) directory which gives an overview of how the intended usage looks like.

## Entrez credentials
`metadatamapping` retrieves data from the Entrez eUtilities using the [`biopython` interface](https://biopython.org/docs/1.75/api/Bio.Entrez.html). By default the Entrez API only allows 3 requests per second if `Entrez.email` and `Entrez.api_key` are not set. This can be increased by setting these properties accordingly which also speeds up the most timeconsuming part of the pipeline which is the accession -> SRA UID mapping as this relies on eSearch which only allows for one accession at a time (maybe it also takes several but I did not test this as I expect it to be cumbersome to pull apart then). So please make sure to set the `Entrez` properties accordingly like so
```python
from Bio import Entrez
Entrez.email = "<user>@<provider>.<domain>"
Entrez.api_key = "<NCBI API key>
```
The email typically is the email associated to your NCBI account. The API key can be generated as described [here](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us)
