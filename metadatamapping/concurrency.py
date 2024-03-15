import os

import multiprocessing as mp
import itertools as it
import pandas as pd

from functools import partial
from typing import Iterable, Callable, Any, Union
from os import PathLike


def read_header(filename: Union[PathLike, str], sep: str) -> list[str]:
    with open(filename, 'r') as file:
        header = file.readline().rstrip().split(sep)
    
    return header


def sequential_writer(
    table: pd.DataFrame, 
    outfilename: Union[PathLike, str]
) -> None:
    header = False if os.path.exists(outfilename) else True

    # this ensures that columns are well aligned between writes
    if header:
        columns = read_header(outfilename, sep = '\t')
        template = pd.DataFrame(columns = columns)
        table = pd.concat([template, table])

    table.to_csv(
        outfilename,
        sep = '\t',
        mode = 'a',
        index = False,
        header = header
    )


def concurrent_writer(
    table: pd.DataFrame, 
    outfilename: Union[PathLike, str], 
    filelock: mp.Manager().Lock
) -> None:
    with filelock:
        sequential_writer(table, outfilename)


def multiprocess_map(
    chunks: Iterable, 
    func: Callable, 
    n_processes: int, 
    function_writes_file: bool,
    **kwargs
) -> list[Any]:
    """
    takes an iterable containing chunks of a whole and a function to process these chunks and uses n_processes to do it
    
    :param chunks:                  Iterable containing chunks of a whole that need to be processed concurrently
    :param func:                    function to process the chunks with
    :param n_processes:             number of concurrent processes to use for processing
    :param function_writes_file:    
    :param **kwargs:                any keyword arguments to pass to func
    
    :return:                        list of anything that is returned by `func`
    """
    with mp.Pool(n_processes) as p:
        if function_writes_file:
            # lock needs to be a managed one otherwise passing it to the pool will fail
            filelock = mp.Manager().Lock()
            map_function = partial(
                func, 
                filelock = filelock, 
                **kwargs
            )
        
        else:
            map_function = partial(
                func,
                **kwargs
            )

        # this needs to be invoked by iterating over it
        results = [result for result in p.imap(map_function, chunks)]
    
    return results
    
            
def singleprocess_map(chunks: Iterable, func: Callable, **kwargs) -> list[Any]:
    """
    processes chunks with func using map. single threaded equivalent to multiprocess_map
    
    :param chunks:        Iterable containing chunks of a whole that need to be processed concurrently
    :param func:          function to process the chunks with
    :param **kwargs:      any keyword arguments to pass to func
    
    :return:              list of anything that is returned by `func`
    """
    map_function = partial(func, **kwargs)
    
    results_iterable = map(map_function, chunks)
    # this needs to be invoked by iterating over it
    return list(results_iterable)


def process_data_in_chunks(
    iterable: Iterable, 
    func: Callable, 
    chunksize: int = 5000, 
    n_processes: int = 1, 
    function_writes_file: bool = True,
    **kwargs
) -> list[Any]:
    """
    uses func to process the given iterable in chunks of size chunksize possibly concurrently
    
    :param iterable:                iterable containing the data to be processed with func
    :param func:                    function that processes the contents of iterable
    :param chunksize:               size of the indiviually processed chunks of the input iterable
    :param n_processes:             number of processes to use for mapping if n_processes > 1 this will be done concurrently using multiprocessing.imap
    :param function_writes_file:    if True multiplrocessing.Manager.Lock is passed to the function with the filelock keyword
                                    (make sure to set this to False if your function does not write a file and you use the concurrent processing)
    :param **kwargs:                any keyword arguments that need to be passed to func
    
    :return:                        list of anything that is returned by `func`
    """
    chunks = it.batched(
        iterable,
        n = chunksize
    )

    if n_processes > 1:
        results = multiprocess_map(
            chunks,
            func,
            n_processes,
            function_writes_file,
            **kwargs
        )

    else:
        results = singleprocess_map(
            chunks,
            func,
            **kwargs
        )

    return results