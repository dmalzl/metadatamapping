class ResponseNotOKError(Exception):
    """
    custom exception to reflect something went wrong with
    retrieving metadata from GEO FTP
    """
    def __init__(self, message):
        super().__init__(message)