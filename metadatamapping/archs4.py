import h5py
import gc
import scipy

import pandas as pd
import numpy as np
import anndata as ad

from os import PathLike
from . import concurrency
from . import sparse
from typing import Union, Iterable, Any


def parse_table(archs4_file: Union[str, PathLike], root_key: str, table_key: str) -> dict[str, Any]:
    """
    parses the archs4 h5 file and extracts the data at /root_key/table_key into a dictionary
    see https://maayanlab.cloud/archs4/help.html for more infos on the available options of
    root_key and table_key
    
    :param archs4_file:   path to the archs4 h5 file to parse
    :param root_key:      root key to use for parsing the file metadata
    :param table_key:     table key to use for parsing the file metadata
    
    :return:              dictionary containing the metadata with column names as keys and column contents as values
    """
    archs4 = h5py.File(archs4_file, 'r')
    
    data = dict()
    table = archs4[root_key][table_key]
    for key in table.keys():
        column_data_piece = table[key][0]
        if isinstance(column_data_piece, np.number):
            column_data = table[key][:]
        
        else:
            column_data = [x.decode('utf-8') for x in table[key][:]]
        
        data[key] = column_data
    
    return data


def filter_table_data(data: dict, retain_keys: Iterable[str]) -> dict[str, Any]:
    """
    filters the data dictionary to retain only the keys specified by retain_keys
    
    :param data:          dictionary as returned by parse_table
    :param retain_keys:   list of data keys to retain
    
    :return:              filtered data dictionary
    """
    retain_data = {
        key: data[key].copy() for key in retain_keys
    }
    return retain_data


def get_srx_samn_from_items(relation_items: Iterable[str]) -> tuple[str, str]:
    """
    takes list of strings generated from an item in the 'relation' column
    of the /meta/samples table in the archs4 h5 file and extracts the SRA and
    BioSample accessions from it
    
    :param relation_items:   list of strings of the format 'key: value'
    
    :return:                 sra and bisample accession as strings
    """
    srx, samn = None, None
    for item in relation_items:
        if not item:
            continue
            
        key, value = item.split(': ')
            
        if key == 'SRA':
            srx = value.split('=')[-1]
        
        if key == 'BioSample':
            samn = value.split('/')[-1]
        
    return srx, samn


def extract_srx_and_samn_accessions(relation_data: Iterable[str]) -> dict[str, list[str]]:
    """
    takes the full 'relation' column of the archs4 h5 /meta/samples table and extracts
    SRA and BioSample accessions for each of the items contained
    
    :param relation_data:  list of strings of the format 'key1: value1, key2: value2, ...'
    
    :return:               dictionary containing keys srx_accession, biosample_accession with lists of these accessions as values
    """
    srx_samn_accessions = {
        key: [] for key in ['srx_accession', 'biosample_accession']
    }
    for relation in relation_data:
        relation_items = relation.split(',')
        srx, samn = get_srx_samn_from_items(
            relation_items
        )
        srx_samn_accessions['srx_accession'].append(srx)
        srx_samn_accessions['biosample_accession'].append(samn)
    
    return srx_samn_accessions


def get_filtered_sample_metadata(archs4_file: Union[str, PathLike], keys_to_retain: Iterable[str]) -> pd.DataFrame:
    """
    parses the archs4 h5 file and returns a pandas.DataFrame of with all columns of the /meta/samples table 
    specified in keys_to_retain plus extracted SRA and BioSamples accessions (columns srx_accession and biosample_accession)
    
    :param archs4_file:       path to the archs4 file to parse
    :param keys_to_retain:    list of keys in the /meta/samples table to retain
    
    :return:                  pandas.DataFrame containing the respective metadata + SRA and BioSample accessions
    """
    data = parse_table(
        archs4_file,
        'meta',
        'samples'
    )

    filtered_data = filter_table_data(
        data, 
        keys_to_retain
    )

    srx_samn_accessions = extract_srx_and_samn_accessions(
        filtered_data['relation']
    )

    srx_samn_table = pd.DataFrame().from_dict(
        srx_samn_accessions,
        orient = 'columns'
    )

    table = pd.DataFrame().from_dict(
        filtered_data, 
        orient = 'columns'
    )

    table = pd.concat(
        [table, srx_samn_table],
        axis = 1
    )

    del data, filtered_data, srx_samn_accessions
    gc.collect()
    
    return table


# adapted from archs4py as this would not install due to issues with Python 3.12 and numpy requirement
def load_data(
    file: Union[PathLike, str], 
    sample_idx: Iterable[int], 
    gene_idx: Iterable[int] = [], 
    n_processes: int = 1
) -> ad.AnnData:
    """
    Retrieve gene expression data from a specified file for the given sample and gene indices.

    :param file:            the file path or object containing the data.
    :param sample_idx:      a list of sample indices to retrieve expression data for.
    :param gene_idx:        a list of gene indices to retrieve expression data for. Defaults to an empty list (return all).

    :return:                AnnData object containing the expression data
    """
    sample_idx = sorted(sample_idx)
    gene_idx = sorted(gene_idx)
    row_encoding = get_gene_name_column(file)
    with h5py.File(file, "r") as f:
        genes = np.array([x.decode("UTF-8") for x in np.array(f[row_encoding])])
        if len(sample_idx) == 0:
            return pd.DataFrame(index=genes[gene_idx])
        gsm_ids = np.array([x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])])[sample_idx]

    if len(gene_idx) == 0:
        gene_idx = list(range(len(genes)))

    data = concurrency.process_data_in_chunks(
        zip(sample_idx, gsm_ids),
        load_samples_in_consecutive_blocks,
        n_processes = n_processes,
        file = file,
        gene_idx = gene_idx
    )

    # scipy.sparse.vstack uses a lot of memory when concatenating the data chunks
    # estimated to double/triple the size of the originally loaded data
    # therefore we use an optimized version of it that is only supposed to concatenate
    # csr_matrices
    sparse_exp = sparse.bmat(data)
    sparse_exp.eliminate_zeros()

    del data
    gc.collect()

    exp = ad.AnnData(
        X = sparse_exp,
        var = pd.DataFrame(index = genes[gene_idx]),
        obs = pd.DataFrame(index = gsm_ids)
    )
    exp.var_names_make_unique()
    return exp


def consecutive(
    sequence_array: np.ndarray[int], 
    *additional_arrays_to_split: np.ndarray[Any]
) -> list[tuple(np.ndarray)]:
    """
    finds all sequences of consecutive elements in a numpy array
    and splits the input arrays accordingly

    :param sequence_array:                  numpy.ndarray containig a sequence of integers
    :param *additional_arrays_to_split:     any number of additional arrays to split into same chunks as sequence array

    :return:                                list of tuples of chunks of sequence_array and additional arrays            
    """
    split_idx = np.where(np.diff(sequence_array) != 1)[0]+1
    arrays = [sequence_array, *additional_arrays_to_split]
    return [np.split(array, split_idx) for array in arrays]


def get_gene_name_column(file: Union[PathLike, str]) -> str:
    """
    detect the conatined gene/transcript names and return corresponding h5 path

    :param file:        string denoting the path to the ARCHS4 file

    :return:            string giving the h5 path to gene/transcript names
    """
    with h5py.File(file) as f:
        if "genes" in list(f["meta"].keys()):
            if "gene_symbol" in list(f["meta/genes"].keys()):
                return "meta/genes/gene_symbol"
            elif "symbol" in list(f["meta/genes"].keys()):
                return "meta/genes/symbol"
        elif "transcripts" in list(f["meta"].keys()):
            if "ensembl_id" in list(f["meta/trancripts"].keys()):
                return "meta/trancripts/ensembl_id"
        else:
            raise Exception("error in gene/transcript meta data")


def read_sample_data(
    file: Union[PathLike, str], 
    sample_idx: Iterable[int], 
    gene_idx: Iterable[int]
) -> scipy.sparse.csr_matrix:
    """
    reads a given block of consecutive samples from the ARCHS4 file

    :param file:                string denoting the path to the ARCHS4 file
    :param sample_idx:          iterable of integers denoting the samples to load from file
    :param gene_idx:            iterable of integers containing the indexes of genes to retain

    :return:                    scipy.sparse.csr_matrix containing the loaded expression data
    """
    with h5py.File(file, "r") as f:
        dense_expression = f["data/expression"][:, sample_idx][gene_idx, :]
        sparse_expression = scipy.sparse.csr_matrix(dense_expression.T)
        sparse_expression.eliminate_zeros()

        del dense_expression
        gc.collect()

    return sparse_expression
    

# this is a wrapper to ensure compliance with API
def load_samples_in_consecutive_blocks(
    sample_idx_gsms: Iterable[tuple[int, str]], 
    file: Union[PathLike, str], 
    gene_idx: Iterable[int], 
    **kwargs
) -> ad.AnnData:
    """
    load the data in consecutive blocks of samples to minimize IO

    :param sample_idx_gsms:     iterable of tuples containing the index of the sample and its accession
    :param file:                string denoting the path to the ARCHS4 file
    :param gene_idx:            iterable of integers containing the indexes of genes to retain
    :param **kwargs:            Any other keyword arguments. This does not have an effect but is just to comply to
                                'process_data_in_chunks' API

    :return:                    AnnData object containing the data
    """
    sample_idxs = np.array([item[0] for item in sample_idx_gsms])
    sample_gsms = np.array([item[1] for item in sample_idx_gsms])
        
    data = [
        read_sample_data(file, consecutive_idx, gene_idx) 
        for consecutive_idx, _ in zip(*consecutive(sample_idxs, sample_gsms))
    ]
    return sparse.bmat(data)


def samples(file, sample_ids, silent=False, n_processes = 1):
    sample_ids = set(sample_ids)
    with h5py.File(file, "r") as f:
        samples = [x.decode("UTF-8") for x in np.array(f["meta/samples/geo_accession"])]

    idx = [i for i,x in enumerate(samples) if x in sample_ids]
    if len(idx) > 0:
        return load_data(file, idx, n_processes=n_processes)
    