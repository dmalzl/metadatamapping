import requests
import logging 

import pandas as pd

from . import dbutils
from . import concurrency
from . import parsers


logging.basicConfig(
    format = '%(threadName)s: %(asctime)s-%(levelname)s-%(message)s',
    datefmt = '%Y-%m-%d %H:%M:%S',
    level = logging.INFO
)

# need to fetch from ftp because e utilities do not provide functionality
# can be done directly from archs4 using gsm and gse accessions
def fetch_soft_metadata(accession: str) -> list[str]:
    """
    fetches metadata from the GEO FTP server. Returns a (possibly) empty list

    :param accession:   string containing the GEO accession to fetch metadata for

    :return:            list of SOFT formatted strings if retrieval was successful else empty list
    """
    URL_BY_ACC = {
        "GSE": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=gse&acc={ACCESSION}&form=text&view=full",
        "GSM": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=gsm&acc={ACCESSION}&form=text&view=full",
    }
    url_base = URL_BY_ACC[accession[:3]]
    response = requests.get(
        url_base.format(ACCESSION = accession)
    )
    
    if not response.ok:
        logging.info(f'could not retrieve {accession}')
        return []

    response.encoding = "UTF-8"
    result_text = response.text
    result_list = result_text.replace("\r", "").split("\n")
    result_list = [elem for elem in result_list if len(elem) > 0]
        
    return result_list


def fetch_and_parse_gse_softs(
    gse_accessions: list[str], 
    already_fetched_softs: dict[str, list[str]]
) -> dict[str, str]:
    """
    fetches GEO Series metadata from GEO FTP. only returns the metadata of the
    oldest Series entry to make sure we retreive the original description.

    :param gse_accessions:          list of GEO series accessions as strings
    :param already_fetched_softs:   dictionary of already fetched GSE metadata to avoid unnecessary API calls

    :return:                        dictionary containing the metadata of interest of the oldest GEO series retrieved
    """
    gse_keys = {
        'Series_summary': 'series_summary',
        'Series_overall_design': 'series_design',
        'Series_submission_date': 'series_submission_date'
    }
    original_series_metadata = {}
    oldest_submission_date = None
    for gse in gse_accessions:
        if gse in already_fetched_softs:
            gse_soft = already_fetched_softs[gse]
            
        else:
            gse_soft = dbutils.retry(
                fetch_soft_metadata,
                gse
            )
            already_fetched_softs[gse] = gse_soft
        
        if parsers.is_superseries(gse_soft):
            continue
        
        gse_metadata = parsers.parse_soft_metadata(gse_soft, gse_keys)
        this_submission_date = parsers.to_date(
            gse_metadata['series_submission_date']
        )
        gse_metadata['series_submission_date'] = this_submission_date
        
        if not original_series_metadata:
            original_series_metadata = gse_metadata
            oldest_submission_date = this_submission_date
            continue
        
        if this_submission_date < oldest_submission_date:
            original_series_metadata = gse_metadata
            
    for key, value in original_series_metadata.items():
        if isinstance(value, list):
            original_series_metadata[key] = ''.join(value)
    
    return original_series_metadata


def fetch_geo_metadata(geo_accessions: list[tuple[str, str]]) -> pd.DataFrame:
    """
    retrieves the metadata for a list of GEO accessions and 
    returns it as a pandas.DataFrame

    :param geo_accessions:  list of GSM, GSE accession tuples. 
                            GSE can be a string of multiple GSE accessions separated by a ','

    :return:                pandas.DataFrame of retrieved metadata
    """
    gsm_keys = {
        'Sample_treatment_protocol_ch1': 'treatment_protocol',
        'Sample_growth_protocol_ch1': 'growth_protocol'
    }

    metadata = {}
    already_fetched_gses = {}
    for gsm, gse in geo_accessions:
        gsm_soft = dbutils.retry(
            fetch_soft_metadata,
            gsm
        )
        gsm_metadata = parsers.parse_soft_metadata(gsm_soft, gsm_keys)
        gsm_metadata.update(
            fetch_and_parse_gse_softs(
                gse.split(','),
                already_fetched_gses
            )
        )
        metadata[gsm] = gsm_metadata

    metadata = pd.DataFrame.from_dict(
        metadata,
        orient = 'index'
    )

    return metadata


def geo_metadata(
    geo_accessions: list[tuple[str, str]], 
    n_processes: int = 1,
    chunksize = 10000
) -> pd.DataFrame:
    """
    retrieves the metadata for a list of GEO accessions and 
    returns it as a pandas.DataFrame. More or less a concurrent version
    of `fetch_geo_metadata`

    :param geo_accessions:  list of GSM, GSE accession tuples. 
                            GSE can be a string of multiple GSE accessions separated by a ','
    :param n_processes:     how many concurrent processes to use to retrieve the data
    :param chunksize:       size of the chunks that are used for each process.

    :return:                pandas.DataFrame of retrieved metadata
    """
    
    metadata_frames = concurrency.process_data_in_chunks(
        geo_accessions,
        fetch_geo_metadata,
        n_processes = n_processes,
        chunksize = chunksize
    )

    return pd.concat(metadata_frames)
