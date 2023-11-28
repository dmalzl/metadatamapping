import pandas as pd


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
    ref = normalize_string(x.iloc[0])
    return all(normalize_string(item) == ref for _, item in x.items())