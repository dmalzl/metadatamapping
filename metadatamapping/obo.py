import obonet as on


def make_id_to_term_dictionary(obo_file_name: str) -> dict[str, str]:
    """
    uses obonet to parse the given obo file and returns a dictionary
    containing ontology ID to onotology term mappings

    :param obo_file_name:   string containing the path to the obo file to read

    :return:                dictionary containing ID to term mapping with ID as key and term as values
    """

    obo_graph = on.read_obo(obo_file_name)
    return {id_: data.get('name') for id_, data in obo_graph.nodes(data=True)}


def map_id_to_term(ont_id: str, ontologies: list[dict]) -> str:
    """
    maps a given ontology ID to it's human readable term. Returns an empty string if ID is not
    contained in any of the given ID to term mappings

    :param ont_id:      string containing the ontology ID to be mapped to its term
    :param ontologies:  list of id_to_term dictionaries for all used ontologies

    :return:            ontology term or empty string if ID is not found in any of the ontologies
    """
    term = None
    for ontology in ontologies:
        if not ont_id in ontology:
            continue

        term = ontology[ont_id]
    
    return term
