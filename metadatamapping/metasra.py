import certifi
import urllib3
import json

import pandas as pd
import itertools as it

from typing import Iterable, Union
from os import PathLike

def study_id_to_metasra(study_ids: Iterable[str]) -> pd.DataFrame:
    """
    used the MetaSRA API to retrieve normalized metdata for all samples in 
    a study given by study_ids from their database
    
    :param study_ids:   iterable containing SRA Biostudy accessions

    :return:            pandas.DataFrame containing MetaSRA normalised metadata
    """
    metasra_data = []
    for i, study_id in enumerate(study_ids):
        # this was necessary before metasra switched to https
        # http = urllib3.PoolManager(
        #     ca_certs = certifi.where(),
        #     cert_reqs = 'CERT_NONE'                      
        # )

        # retrieving number of studies that match search criteria
        r = urllib3.request(
            'GET', 
            'https://metasra.biostat.wisc.edu/api/v01/samples.csv',
            fields = {
                'study': study_id,
                'species': 'human', 
                'assay': 'RNA-seq', 
                'limit': 10000,
                'skip': 0
            }
        )

        study_data = pd.read_csv(BytesIO(r.data))
        
        if not study_data.empty:
            metasra_data.append(study_data)
    
    return pd.concat(metasra_data)


def convert_attributes_to_dict(attributes: str) -> dict[str, str]:
    """
    converts the attributes string into an attributes dictionary

    :param attributes:  string of key-value pairs of the format key1: value1; key2: value2; ...

    :return:            dictionary of attributes with key as key and value as value
    """
    attribute_dict = dict()
    for kv_pair in attributes.split('; '):
        k, v = kv_pair.split(': ')
        attribute_dict[k] = v
    
    return attribute_dict


def convert_to_dict_and_write_json(metadata_rows: Iterable[pd.Series], output_json: Union[PathLike, str]) -> None:
    """
    takes an iterable of pandas.DataFrame rows as pd.Series with keys 'accession' and 'attribute' converts them into
    a dictionary with key accession and attribute dictionary as values and writes this as JSON to output_json

    :param metadata_rows:   iterable of pd.Series with keys 'accession' and 'attribute'
    :param output_json:     name of the outputfile to write the JSON formatted data to
    
    :return:                None
    """
    json_dict = {
        row.accession: convert_attributes_to_dict(row.attribute) for _, row in metadata_rows
    }
    with open(output_json, 'w') as jsonfile:
        json.dump(json_dict, jsonfile, indent=4, separators=(',', ': '))


def raw_metadata_to_json(metadata: pd.DataFrame, output_json: Union[PathLike, str], chunksize = -1) -> None:
    """
    takes a pandas.DataFrame with columns 'accession' and 'attribute' and converts it to a
    MetaSRA input JSON style dictionary (accession as key, dictionary of attributes as value)
    which is then written to a JSON file. If chunksize > 0 the DataFrame is split into chunks
    and multiple outputfiles will be created.

    :param metadata:        pandas.DataFrame with columns 'accession' and 'attribute'
    :param output_json:     name of the output JSON file, will be appended with numbers 1 - n_chunks if chunksize > 0
    :param chunksize:       number of individual samples per chunk, n_chunks = n_samples / chunksize

    :return:                None
    """
    output_file_prefix = '.'.join(output_json.split('.')[:-1])
    output_file_suffix = output_json.split('.')[-1]
    if chunksize > 0:
        chunks = it.batched(
            metadata.iterrows(),
            chunksize
        )
        for i, chunk in enumerate(chunks):
            convert_to_dict_and_write_json(
                chunk, 
                f'{output_file_prefix}_{i+1}.{output_file_suffix}'
            )

    else:
        convert_to_dict_and_write_json(
            metadata.iterrows(),
            f'{output_file_prefix}.{output_file_suffix}'
        )


def metasra_output_json_to_dataframe(output_json: Union[PathLike, str]) -> pd.DataFrame:
    """
    """
    import json
    with open(output_json, 'r') as f:
        mappings = json.load(f)

    metasra_data = []
    for accession, mapping in mappings.items():
        sample_type = mapping['sample type']
        sample_type_confidence = mapping['sample-type confidence']
        mapped_ontology_terms, mapped_ontology_ids = [], []
        for term in mapping['mapped ontology terms']:
            term_id, term_name = term.split('|')
            mapped_ontology_terms.append(term_name)
            mapped_ontology_ids.append(term_id)
        
        metasra_data.append(
            [
                accession, 
                sample_type, 
                sample_type_confidence, 
                ', '.join(mapped_ontology_ids), 
                ', '.join(mapped_ontology_terms)
            ]
        )

    metasra_df = pd.DataFrame(
        metasra_data, 
        columns = [
            'accession', 
            'sample_type', 
            'sample_type_confidence', 
            'mapped_ontology_ids', 
            'mapped_ontology_terms'
        ]
    )
    return metasra_df
        

def make_accession_and_attributes_table(archs4_annotated: pd.DataFrame, accesseion_column: str = 'sample') -> pd.DataFrame:
    """
    takes an annotated ARCHS4 metdata table as returned by metadata.merge_to_annotated_metadata
    and extracts a given accession column and the raw BioSample attributes. removes samples
    that do not have any attributes

    :param archs4_annotated:    annotated ARCHS4 metadata table (see metadata.merge_to_annotated_metadata)
    :param accession_column:    string denoting the column with the desired NCBI accession

    :return:                    pandas.DataFrame with columns 'accession' and 'attribute'
    """
    accession_and_attributes = archs4_annotated.loc[
        :, 
        [accession_column, 'raw_biosample_metadata']
    ].rename(
        columns = {
            accession_column: 'accession',
            'raw_biosample_metadata': 'attribute'
        }
    ).dropna()

    return accession_and_attributes 
