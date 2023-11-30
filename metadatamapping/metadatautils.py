import pandas as pd
from unidecode import unidecode


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


def normalize_string(string: str, extra_chars_to_remove = ' ') -> str:
    """
    removes non ascii characters from a given string

    :param string:      string from which to remove non ascii characters

    :return:            string with just ascii characters
    """
    normalized_string = []
    for char in string:
        if ord(char) < 128 and not char in extra_chars_to_remove:
            normalized_string.append(char)
    
    return ''.join(normalized_string)


def all_equal(x: pd.Series) -> bool:
    """
    checks if all items in a pandas.Series are the same

    :param x:       pandas.Series object to check

    :return:        True if all items are the same else False
    """
    ref = unidecode(remove_whitespace(x.iloc[0]))
    return all(unidecode(remove_whitespace(item)) == ref for _, item in x.items())
