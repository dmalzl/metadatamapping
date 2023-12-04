import regex

import pandas as pd

from unidecode import unidecode
from typing import Iterable


def remove_whitespace(string: str) -> str:
    """
    removes whitespace from string

    :param string:      string object from which to remove whitespaces

    :return:            string without whitespace
    """
    no_whitespace_string = []
    for char in string:
        if char != ' ':
            no_whitespace_string.append(char)
    
    return ''.join(no_whitespace_string)


def all_equal(x: pd.Series) -> bool:
    """
    checks if all items in a pandas.Series are the same

    :param x:       pandas.Series object to check

    :return:        True if all items are the same else False
    """
    ref = unidecode(remove_whitespace(x.iloc[0]))
    return all(unidecode(remove_whitespace(item)) == ref for _, item in x.items())


def remove_empty_string_tokens(tokens: list[str]) -> list[str]:
    return [token for token in tokens if token]


def normalize_term(term: str, first_token_replacement: str) -> str:
    """
    this function fixes spelling mistakes and relies on the 
    assumption that term contains the word to be replaces at 
    the first position of the word sequence. Additionally, also
    replaces '-', '/' and '_' with ' '

    :param term:                        string containing the term
    :param first_token_replacement:     string used as a replacement for the first word in term

    :return:                            string containing the normalize term
    """
    replace_strings = ['-', '/', '_']
    for replace_string in replace_strings:
        term = term.replace(replace_string, ' ')
        
    tokens = remove_empty_string_tokens(term.split(' '))
    tokens[0] = first_token_replacement
    return ' '.join(tokens)


def extract_treatment(raw_metadata_items: list[str], matchers: dict[regex.Pattern]) -> str:
    """
    extracts all normalized terms and their values that match any of the regex patterns in matchers

    :param raw_metadata_items:  list of strings of the form "term: value"
    :param matchers:            dictionary containing first_token_replacement as key and regex.Patterns for matching as values

    :return:                    string of matched "term: value" pairs concatenated with '; ' or empty string if nothing matches
    """
    extracted_items = []
    for raw_metadata_item in raw_metadata_items:
        for term_name, term_matcher in matchers.items():
            m = term_matcher.match(raw_metadata_item)
            if not m:
                continue

            term, value = raw_metadata_item.split(': ')
            term = normalize_term(term, term_name)
            extracted_items.append(f'{term}: {value}')

            # guards against possible double matches 
            # which technically should not occur
            break
    
    return '; '.join(extracted_items)


def map_treatment(raw_sample_metadata: Iterable[str]) -> list[str]:
    # the regular expressions in this dictionary were built by non-exhaustive
    # trial and error and thus might not match everything that is contained in
    # the raw metadata but should match most of it. This trade-off is made
    # with the hope that the shared embedding of the model that will be trained on this
    # will anyway bring these samples together
    term_matchers = {
        'treatment': regex.compile('^t(reatm){i<=1,d<=1,s<=1,2i+2d+1s<=4}'),
        'transfection': regex.compile('^t(ransfect){i<=1,d<=1,s<=1,2i+2d+1s<=4}'),
        'transduction': regex.compile('^t(ransduc){i<=1,d<=1,s<=1,2i+2d+1s<=4}')
    }

    mapped_treatments = []
    for raw_metadata in raw_sample_metadata:
        raw_metadata_items = raw_metadata.split('; ')
        mapped_treatments.append(
            extract_treatment(raw_metadata_items, term_matchers)
        )
            
    return mapped_treatments
