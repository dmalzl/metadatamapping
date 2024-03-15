import logging

import pandas as pd

from . import parsers
from . import summary
from . import concurrency
from . import dbutils
from . import metasra
from . import metadatautils
from . import geo
from typing import Union, Iterable
from os import PathLike

logging.basicConfig(
    format = '%(threadName)s: %(asctime)s-%(levelname)s-%(message)s',
    datefmt = '%Y-%m-%d %H:%M:%S',
    level = logging.INFO
)

def map_accessions_to_srauids(
    table: pd.DataFrame, 
    outfilename: Union[PathLike, str], 
    chunksize: int = 5000, 
    n_processes: int = 1
) -> None:
    """
    takes a pandas.DataFrame containing columns srx_accession, biosample_accession and geo_accession, maps those accessions
    to SRA UIDs and writes the resulting map to outfilename. If n_processes > 1 this will be done concurrently
    
    :param table:         pandas.DataFrame containing columns srx_accession, biosample_accession and geo_accession
    :param db:            string denoting the NCBI database to retrieve the UIDs from
    :param outfilename:   path to the outputfile the resulting map should be written to
    :param chunksize:     size of the indiviually processed chunks of the input table
    :param n_processes:   number of processes to use for mapping if n_processes > 1 this will be done concurrently using multiprocessing.imap
    
    :return:              None
    """
    mapping_table = parsers.get_not_yet_mapped_uids(
        table,
        outfilename,
        'accession',
        ['geo_accession', 'srx_accession', 'biosample_accession']
    )
    
    concurrency.process_data_in_chunks(
        mapping_table.iterrows(),
        summary.map_accessions_to_uids,
        db = 'sra',
        outfilename = outfilename,
        n_processes = n_processes,
        chunksize = chunksize
    )


def determine_accession_type(accession: str) -> str:
    """
    determines the type of accession used for retrieving the UIDs
    currently identifies SRX SAMN or GSM accessions

    :param accession:   string containing the accession to be identified

    :return:            string giving the type of accession
    """
    types = {
        'experiment': 'SRX',
        'biosample': 'SAMN',
        'geo': 'GSM'
    }
    
    for accession_type, accession_prefix in types.items():
        if accession.startswith(accession_prefix):
            return accession_type
        

def merge_uids_and_accessions(uids: pd.DataFrame, accessions: pd.DataFrame) -> pd.DataFrame:
    """
    merges the uid frame with the retrieved accessions

    :param uids:            pandas.DataFrame containing 'uid' and 'accession' columns
    :param accessions:      pandas.DataFrame containing accessions correspondig to uids

    :return:                pandas.DataFrame merged on uid
    """
    uids = uids.copy()
    uids['accession_type'] = uids.accession.apply(
        determine_accession_type
    )

    merged_groups = []
    for accession_type, group in uids.groupby('accession_type'):
        group = group.rename(
            columns = {'accession': accession_type}
        )
        
        merged_group = group.merge(
            accessions,
            on = accession_type,
            how = 'inner'
        )

        merged_groups.append(merged_group)
        
    return pd.concat(merged_groups)


def srauids_to_accessions(srauids: pd.DataFrame, chunksize: int = 50000) -> pd.DataFrame:
    """
    retrieves all available accessions for the sra uids given in 'uid' column of srauids

    :param srauids:     pandas.DataFrame containing an 'accession' and a 'uid' column
    :param chunksize:   integer giving the size of the chunks posted and retrieved from Entrez

    :return:            pandas.DataFrame containing experiment, sample, biosample, study, bioproject and geo accessions
    """
    accessions = summary.summaries_from_uids(
        srauids.uid,
        'sra',
        parsers.accessions_from_esummary_response,
        parsers.accession_matchers,
        chunksize = chunksize
    )
    
    return merge_uids_and_accessions(srauids, accessions)


def biosample_uids_to_metadata(biosample_uids: Iterable[Union[str, int]], chunksize: int = 50000) -> pd.DataFrame:
    """
    retrieves all available metadata for the biosample uids

    :param biosample_uids:  iterable of biosample UIDs as string or int
    :param chunksize:       integer giving the size of the chunks posted and retrieved from Entrez

    :return:            pandas.DataFrame containing experiment, sample, biosample, study, bioproject and geo accessions
    """
    metadata = summary.summaries_from_uids(
        biosample_uids,
        'biosample',
        parsers.metadata_from_biosample_uids,
        parsers.biosample_metadata_parsers,
        chunksize = chunksize
    )
    
    return metadata


def metasra_from_study_id(study_ids: Iterable[str]) -> pd.DataFrame:
    """
    used the MetaSRA API to retrieve normalized metdata for all samples in 
    a study given by study_ids from their database
    
    :param study_ids:   iterable containing SRA Biostudy accessions

    :return:            pandas.DataFrame containing MetaSRA normalised metadata
    """
    metasra_data = dbutils.process_in_chunks(
        study_ids,
        metasra.study_id_to_metasra,
        chunksize = 1000
    )

    return pd.concat(metasra_data)


def merge_to_annotated_metadata_frame(
    biosample_metadata: pd.DataFrame, 
    srauids_to_biosampleuids: pd.DataFrame,
    ncbi_accessions: pd.DataFrame,
    archs4_metadata: pd.DataFrame,
    geo_metadata: pd.DataFrame
) -> pd.DataFrame:
    """
    combines the passed pandas.DataFrames to a fully annotated ARCHS4 metadata set

    :param biosample_metadata:          pandas.DataFrame containing BioSample metadata (see biosample_uids_to_metadata)
    :param srauids_to_biosampleuids:    pandas.DataFrame containing the linked SRA to BioSample UIDs (see link.link_sra_to_biosample)
    :param ncbi_accessions:             pandas.DataFrame containing the SRA accessions (see srauids_to_accessions)
    :param archs4_metadata:             pandas.DataFrame containing the extracted ARCHS4 metadata (see archs4.get_filtered_sample_metadata)
    :param geo_metadata:                pandas.DataFrame containing the fetched GEO metadata (see geo.geo_metadata)

    :return:                            pandas.DataFrame containing annotated ARCHS4 metadata
    """
    biosample_uids_and_metadata = srauids_to_biosampleuids.merge(
        biosample_metadata,
        on = 'biosample',
        how = 'inner'
    )

    accessions_biosample_metadata_uids = biosample_uids_and_metadata \
        .rename(
            columns = {
                'sra': 'sra_uid',
                'biosample': 'biosample_uid',
                'title': 'biosample_title',
                'attribute': 'raw_biosample_metadata'
            }
        ).merge(
            ncbi_accessions.rename(
                columns = {
                    'uid': 'sra_uid'
                }
            ),
            on = 'sra_uid',
            how = 'inner'
        )

    archs4_annotated = accessions_biosample_metadata_uids.merge(
        archs4_metadata.drop(
            columns = [
                    'srx_accession', 
                    'biosample_accession'
                ]
        ).rename(
            columns = {
                'geo_accession': 'geo',
                'title': 'geo_title',
                'source_name_ch1': 'geo_source_name',
                'characteristics_ch1': 'geo_metadata'
            }
        ),
        on = 'geo',
        how = 'inner'
    )

    # this is done in order to remove scrambled annotations due to spurious links to BioSample
    spurious_duplication = archs4_annotated.sra_uid.duplicated(keep = False)
    archs4_annotated_nonduplicated = archs4_annotated.loc[~spurious_duplication, :]
    archs4_annotated_duplicated = archs4_annotated.loc[spurious_duplication, :]

    all_equal_idx = archs4_annotated_duplicated[['geo_title', 'biosample_title']].apply(
        metadatautils.all_equal,
        axis = 1
    )

    archs4_annotated_deduplicated = pd.concat(
        [
            archs4_annotated_nonduplicated,
            archs4_annotated_duplicated.loc[all_equal_idx, :]
        ]
    )

    archs4_annotated_deduplicated = archs4_annotated_deduplicated.merge(
        geo_metadata,
        on = 'geo',
        how = 'inner'
    )
    return archs4_annotated_deduplicated.reset_index(drop = True)


def annotated_with_metasra_and_treatment(annotated_frame: pd.DataFrame, metasra_frame: pd.DataFrame) -> pd.DataFrame:
    """
    merges the annotated DataFrame returned by 'merge_to_annotated_metadata_frame' with the normalized metadata DataFrame
    produced by MetaSRA and as returned by 'metasra.metasra_output_json_to_dataframe'. Additionally, the function also 
    attempts to parse out any treatment information from the raw metadata which is usually not kept by metasra into a
    new column 'treatment'. This includes any raw metadata property that starts with 'treat', 'transfect' or 'transduc'.
    Note that resulting matches are normalized by replacing the first word in the match with either 'treatment', 'transfection'
    or 'transduction' respectively. This is a crude attempt to remove any spelling mistakes.

    :param annotated_frame:     pandas.DataFrame as returned by 'merge_to_annotated_metadata_frame'
    :param metasra_frame:       pandas.DataFrame as returned by 'metasra.metasra_output_json_to_dataframe'

    :return:                    merged pandas.DataFrame containing an additional new column 'treatment'
    """
    sample_accessions = annotated_frame.loc[:, ['sample', 'geo']]
    metasra_accession_annotated = metasra_frame.merge(
        sample_accessions.rename(columns = {'sample': 'accession'}),
        on = 'accession',
        how = 'inner'
    )
    metasra_accession_annotated.drop_duplicates(inplace = True)

    metasra_annotated = metasra_accession_annotated.merge(
        annotated_frame.rename(
            columns = {
                'geo_accession': 'geo',
                'title': 'geo_title',
                'source_name_ch1': 'geo_source_name',
                'characteristics_ch1': 'geo_metadata'
            }
        ),
        on = 'geo',
        how = 'inner'
    )

    metasra_annotated['treatment'] = metadatautils.map_treatment(
        metasra_annotated.raw_biosample_metadata
    )

    return metasra_annotated


def geo_metadata(
    geo_accessions: pd.DataFrame, 
    outfilename: Union[PathLike, str],
    n_processes: int = 1,
    chunksize = 10000
) -> pd.DataFrame:
    """
    retrieves the metadata for a list of GEO accessions and 
    returns it as a pandas.DataFrame. More or less a concurrent version
    of `fetch_geo_metadata`

    :param geo_accessions:  pandas.DataFrame containing at least two columns GSM and GSE 
                            GSE can be a string of multiple GSE accessions separated by a ','
    :param outfilename:     name of the file the retrieved metadata should be written to
    :param n_processes:     how many concurrent processes to use to retrieve the data
    :param chunksize:       size of the chunks that are used for each process.

    :return:                pandas.DataFrame of retrieved metadata
    """
    mapping_table = parsers.get_not_yet_mapped_accessions(
        geo_accessions,
        outfilename,
        'GSM',
        ['geo_accession']
    )
    concurrency.process_data_in_chunks(
        mapping_table.iterrows(),
        geo.fetch_geo_metadata,
        outfilename = outfilename,
        n_processes = n_processes,
        chunksize = chunksize
    )