#! /usr/bin/env python
__author__ = 'BinglanLi'

import json
import sys

import numpy as np
from pathlib import Path
from typing import List, Set, Tuple, Optional

# todo: convert sys.exit(1) to reportableExceptions

# define the list of wobble genotypes
_wobble_genotype_list: list = ['S', 'Y', 'M', 'K', 'R', 'W', 'V', 'H', 'D', 'B', 'N']
# build a dictionary of wobble genotypes and their corresponding basepairs
_wobble_match_table: dict[str, List[str]] = {
    'M': {'A', 'C'},
    'R': {'A', 'G'},
    'W': {'A', 'T'},
    'S': {'C', 'G'},
    'Y': {'C', 'T'},
    'K': {'G', 'T'},
    'V': {'A', 'C', 'G'},
    'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'},
    'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'}
}


def read_json(file: Path) -> dict[str, str]:
    """
    Read json files
    :param file: path to the input json file
    :return: json_data: a dictionary of json data
    """
    try:
        # open and read json files
        with open(file, 'r') as json_file:
            json_data: dict = json.load(json_file)
        return json_data

    except FileNotFoundError:
        print(f"File '{file}' not found.")
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")


def get_allele_defining_variants(json_data: dict) -> dict[str, dict[str, str]]:
    """
    Extract information of allele-defining positions
    :param json_data: a dictionary of allele-defining variants and named alleles from PharmCAT JSON
    :return: a dictionary of {hgvs_name: {
        position: genomic_coordinate,
        ref: reference_genotype} }
    """
    # extract desirable information of allele-defining positions from the json data
    try:
        hgvs_names: List[str] = [entry['chromosomeHgvsName'] for entry in json_data['variants']]
        positions: List[str] = [entry['position'] for entry in json_data['variants']]
        reference_genotypes: List[str] = [entry['ref'] for entry in json_data['variants']]

        # check whether there are any empty values
        if None in hgvs_names:
            print('One of the HGVS names is empty.')
        if None in positions:
            print('Genomic coordinate is empty for one of the allele definition positions.')
        if None in reference_genotypes:
            print('Reference genotype is empty for one of the allele definition positions.')

        # generate a dictionary of all allele-defining positions and its reference genotypes
        allele_defining_variants: dict[str, dict[str, str]] = {
            h: {
                'position': p,
                'ref': r
            } for h, p, r in zip(hgvs_names, positions, reference_genotypes)}

        return allele_defining_variants

    except Exception as e:
        # print warnings if any of the attributes does not exist in the json file
        print(f'An error occurred in one of the {e}.'
              f'Check whether the attributes exist for all allele-defining positions.')


def get_allele_definitions(json_data: dict,
                           variants: dict[str, dict[str, str]]) -> dict[str, dict[list[str], list[str]]]:
    """
    Extract allele definitions
    :param json_data: a dictionary of allele-defining variants and named alleles from PharmCAT JSON
    :param variants: a dictionary of all allele-defining positions
    :return: a dictionary of {allele_name: {
        defining_genotype: list of defining genotypes, whose sequence should correspond to json - 'variants',
        reference_genotype: list of reference genotype at each position} }
    """
    # extract named allele definitions
    try:
        allele_list: List[str] = [entry['name'] for entry in json_data['namedAlleles']]
        geno_list: List[List[str]] = [entry['alleles'] for entry in json_data['namedAlleles']]

        # check whether there are any empty values
        if None in allele_list:
            print('One of the allele names is missing')
        # check whether any allele definition is empty
        if any(x.count(None) == len(x) for x in geno_list):
            print('One of the alleles is not defined by any genotypes at any positions.')
        # check whether each allele is defined by the same number of positions in the json 'variants' attribute
        n_variants = len(variants)
        if any(len(g) != n_variants for g in geno_list):
            print('One of the alleles has a different number of allele-defining positions.')

        # generate a dictionary of allele definitions with the allele name, defining genotypes, and reference genotypes
        reference_genotypes = [x['ref'] for x in variants.values()]
        n_alleles = len(allele_list)
        allele_definitions: dict[str, dict[list[str], list[str]]] = {
            allele_list[i]: {
                'genotypes': geno_list[i],
                'references': reference_genotypes
            } for i in range(n_alleles)
        }

        return allele_definitions

    except Exception as e:
        # print warnings if any of the attributes does not exist in the json file
        print(f'An error occurred in one of the {e}.'
              f'Check whether the attributes exist for all allele-defining positions.')


def count_allele_scores(allele_definitions: dict[str, dict[list[str], list[str]]]) -> dict[str, int]:
    """
    calculate the haplotype scores by counting the number of non-empty allele-defining positions
    :param allele_definitions: a dictionary of allele definitions {n, {p, p} }
    :return: hap_score: dict[str, int] = { allele_name: allele_score }
    """
    # for each allele, count the number of allele-defining genotypes
    n_defining_genotypes: list[int] = [len(x['genotypes']) - x['genotypes'].count(None)
                                       for x in allele_definitions.values()]

    # create a dictionary of haplotype scores where dict['allele name'] = {haplotype score}
    allele_scores: dict[str, int] = dict(zip(allele_definitions, n_defining_genotypes))

    return allele_scores


def count_diplotype_scores(allele_scores: dict[str, int]) -> dict[str, int]:
    """
    Based on a dictionary of allele/haplotype scores, get all the diplotype combinations and their scores
    :param allele_scores: a dictionary where keys = allele name, values = allele score
    :return: dip_score: a dictionary where keys = the diplotype names and the values = the diplotype scores
    """
    # calculate the diplotype scores
    diplotype_scores: dict[str, int] = {k1 + '/' + k2: v1 + v2
                                        for k1, v1 in allele_scores.items()
                                        for k2, v2 in allele_scores.items()}

    return diplotype_scores


def fill_definitions_with_references(allele_definitions: dict[str, dict[list[str], list[str]]]) \
                                     -> dict[str, dict[list[str], list[str], list[str]]]:
    """
    fill up empty allele-defining positions with the reference genotype
    :return: a list of genotypes where empty/None genotypes have been filled with reference genotypes at the position
    """
    for a, sub_dict in allele_definitions.items():
        # get the list of allele defining genotypes and the reference genotypes at the corresponding positions
        definitions: list[str] = sub_dict['genotypes']
        references: list[str] = sub_dict['references']

        # get the list of definitions where NAs are filled with reference genotypes
        n = len(definitions)
        complete_genotypes: List[str] = [definitions[i] if definitions[i] else references[i] for i in range(n)]
        allele_definitions[a]['complete_genotypes'] = complete_genotypes

    return allele_definitions


def split_wobl_nonwobl_genotypes(genotypes: List[str]) -> Tuple[List[str], List[str]]:
    """
    split wobble and non-wobble genotypes from a list into two separate lists
    :param genotypes: a list of genotypes
    :return: g_w: a list of wobble genotypes
    :return: g_nw: a list of non-wobble genotypes
    """
    # set up empty lists for wobble and non-wobble genotypes
    g_w: List[str] = []
    g_nw: List[str] = []

    # iterate over genotypes to find wobble vs non-wobble genotypes
    for g in genotypes:
        (g_nw, g_w)[g in _wobble_genotype_list].append(g)

    return g_w, g_nw


def replace_wobble(wobl_genotype: str) -> List[str]:
    """
    replace wobble genotypes with basic basepairs A/T/C/G
    :param wobl_genotype: an allele-defining genotype
    :return: a list of genotypes with only A, T, C, G, or indels
    """
    g_flat = _wobble_match_table[wobl_genotype]
    return g_flat


def get_unique_combinations(g1: str, g2: str) -> Set[str]:
    """
    get all possible combinations of genotypes at an allele-defining positions
    :param g1: a list of allele-defining genotypes from allele 1
    :param g2: a list of allele-defining genotypes from allele 2
    :return: a set of unique combinations of genotypes
    """
    # replace wobbles
    g1_ = replace_wobble(g1) if g1 in _wobble_genotype_list else [g1]
    g2_ = replace_wobble(g2) if g2 in _wobble_genotype_list else [g2]

    # get all possible genotype combinations at this positions
    uniq_comb: Set[str] = {'/'.join(sorted([x, y])) for x in g1_ for y in g2_}

    return uniq_comb


def get_diplotype_definition_dictionary(allele_definitions: dict[str, dict[list[str], list[str], list[str]]],
                                        allele_defining_variants: dict[str, dict[str, str]]) \
                                        -> dict[str, dict[str, Set[str]]]:
    """
    obtain a dictionary of diplotype definitions
    :param allele_definitions: a dictionary of allele definitions
    :param allele_defining_variants: a dictionary of allele-defining positions
    :return: a dictionary {diplotype_name: {position: {unique combinations of genotypes at a position} } }
    """

    # initialize variables
    # set up an empty dictionary that stores each diplotype's definition
    diplotype_definitions: dict[str, dict[str, Set[str]]] = {}

    # get the list of alleles
    alleles = [*allele_definitions]
    hgvs_names = [*allele_defining_variants]

    # get the length of the allele list for the loop to iterate through each element
    n_allel = len(allele_definitions)

    # iterate over the allele list
    for i in range(n_allel):
        for j in range(i, n_allel):

            # get the allele names
            a1: str = alleles[i]
            a2: str = alleles[j]

            # specify the diplotype's name, which will be used as the key in dic_mat
            diplotype_name = a1 + '/' + a2

            # check whether the diplotype name is already in the dic_mat dictionary
            if diplotype_name in diplotype_definitions:
                # each diplotype should only be processed once, so each key/diplotype name should be unique
                print(f'Diplotype {diplotype_name} is processed more than once. Something is wrong.')
                sys.exit(1)

            # if no reportable errors, continue to generate the diplotype dictionary
            # first, get the complete defining genotypes without empty values for the alleles
            g1 = allele_definitions[a1]['complete_genotypes']
            g2 = allele_definitions[a2]['complete_genotypes']

            # second, get unique combinations of genotypes at each genetic position
            genotype_combinations = [get_unique_combinations(x, y) for x, y in zip(g1, g2)]

            # third, create a dictionary variable to store the unique genotype combinations at each position
            definition_i: dict[str, Set[str]] = dict(zip(hgvs_names, genotype_combinations))

            # finally, append the diplotype and its genotypes to the dictionary
            diplotype_definitions[diplotype_name] = definition_i
        # end - for j in range(n_allel)
    # end - for i in range(i, n_allel)

    return diplotype_definitions


def convert_dictionary_to_numpy_array(diplotype_definitions: dict[str, dict[str, Set[str]]],
                                      allele_defining_variants: dict[str, dict[str, str]]) -> np.ndarray:
    """
    convert the diplotype definition dictionary to a numpy array
    :param allele_defining_variants: a dictionary of {hgvs_name: {position: [], ref: []}}
    :param diplotype_definitions: {'diplotype name': {'position hgvs name': {unique set of genotypes at a position} } }
    :return: a 3D numpy array [diplotypes, genotype combinations, positions]. here genotype combinations indicates 'A/G'
    """
    # initialize variables
    # set up an empty set to record all the possible genotype combinations across all positions in a gene
    unique_genotypes: set[str] = set()
    # use python Set to find all the unique genotype combinations across positions
    unique_genotypes.update([g for sub_dict in diplotype_definitions.values() for x in sub_dict.values() for g in x])
    # convert the genotype set to a list which has an order
    genotype_list: list[str] = list(unique_genotypes)

    # get the list of diplotype names
    diplotype_names: list[str] = [*diplotype_definitions]

    # get the list of hgvs names for all positions
    hgvs_names: list[str] = [*allele_defining_variants]

    # initialize an empty numpy array with rows and columns
    n_diplotypes: int = len(diplotype_names)
    n_positions: int = len(hgvs_names)
    n_genotypes: int = len(genotype_list)
    diplotype_array: np.ndarray = np.zeros((n_diplotypes, n_genotypes, n_positions), dtype='int8')

    # fill the numpy array with 1 where position(key)-genotype(value) combinations exist
    for diplotype, one_diplotype_definition in diplotype_definitions.items():

        # get the diplotype index for the numpy array
        diplotype_idx = diplotype_names.index(diplotype)

        for position, genotypes in one_diplotype_definition.items():

            # get the position index for the numpy array
            position_idx = hgvs_names.index(position)

            for g in genotypes:
                # get the genotype index for the numpy array
                genotype_idx = genotype_list.index(g)

                # fill the corresponding cell in the numpy array with 1 to denote the definition
                diplotype_array[diplotype_idx, genotype_idx, position_idx] = 1

    return diplotype_array


def is_sharing_definitions(g1: np.ndarray, g2: np.ndarray, axis: Optional[int] = 1) -> np.ndarray[bool]:
    """
    to identify whether g1 (a diplotype's definition) shares definitions with another or any other diplotypes
    :param g1: a numpy array of a diplotype's definition (M, P).
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array for one diplotype (m,p) or multiple diplotypes (D, M, P)
                D = the total number of diplotypes for a gene
    :param axis: the axis at which g1 and g2 will be compared
            if g2 is a single diplotype, use axis = 0
            if g2 is a d x m x p array of multiple diplotypes, use axis = 1
    :return: a numpy array with boolean elements of length (1, ) or (D, )
    """

    # calculate product of g1 and g2
    g_prod = g1 * g2

    # check whether g1 shares any definitions with another diplotype by identifying non-empty columns
    # a non-empty column of g_prod means that g1 and this other diplotype shared a definition at a position
    non_empty_cols = np.any(g_prod != 0, axis=axis)

    # find diplotypes that share definitions with g1 at every allele-defining positions
    status = np.all(non_empty_cols, axis=axis)

    return status


def is_discrepant(g1: np.ndarray, g2: np.ndarray, axis: Optional[int] = 1) -> np.ndarray[bool]:
    """
    to identify whether g1 (a diplotype's definition) is different from another or any other diplotypes
    :param g1: a numpy array of a diplotype's definition (M, P).
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array for one diplotype (m,p) or multiple diplotypes (D, M, P)
                D = the total number of diplotypes for a gene
    :param axis: the axis at which g1 and g2 will be compared
            if g2 is a single diplotype, use axis = 0
            if g2 is a d x m x p array of multiple diplotypes, use axis = 1
    :return: a numpy array with boolean elements of length (1, ) or (D, )
\
    """
    # calculate product of g1 and g2
    g_prod = g1 * g2

    # check whether there is any zero column
    # an empty column means that g1 does not share defining genotypes with a certain diplotype
    # np.all(axis=axis) checks whether a position column has only empty entries
    empty_cols = np.all(g_prod == 0, axis=axis)

    # identify diplotypes that differ from g1 by identifying empty position columns
    status = np.any(empty_cols, axis=axis)

    return status


def is_equivalent(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 and g2 are equivalent
    :param g1: a numpy array of a diplotype's definition (M, P)
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # compare d1 with d2 to find equivalent diplotypes
    # 'np.all(g1 == g2, axis=(1, 2))' is faster than '[numpy.array_equal(g1, g) for g in g2]'
    status = np.all(g1 == g2)

    return status


def is_included(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 is included in g2
    in other words, whether g1 can be considered as a sub-allele of g2

    :param g1: a numpy array of a diplotype's definition (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # use numpy vectorized manipulation to check whether g1 is included as a sub-allele by any other diplotypes
    g_subs = g1 - g2

    # find diplotypes that includes g1
    status = set(g_subs.flatten()) == {0, -1}

    return status


def is_inclusive(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 includes g2
    in other words, whether g1 can be considered as a super-allele of g2

    :param g1: a numpy array of a diplotype's definition (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # use numpy vectorized manipulation to check whether g1 is included as a sub-allele by any other diplotypes
    g_subs = g1 - g2

    # find diplotypes that includes g1
    status = set(g_subs.flatten()) == {0, 1}

    return status


def is_overlapping(g1: np.ndarray, g2: np.ndarray) -> bool:
    """
    to identify whether g1 overlaps with g2
    g1 and g2 must have intersecting definitions, but g1 and g2 cannot be equivalent or a sub-allele of one the other

    :param g1: a numpy array of a diplotype's definition (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: a numpy array of a diplotype's definition (M, P)
    :return: True or False
    """
    # calculate g1*g2 and g1-g2
    g_prod = g1 * g2
    g_subs = g1 - g2

    # make sure g_prod does not have any empty columns
    # an empty column means that g1 and g2 do not share defining genotypes at a position
    non_empty_cols = np.any(g_prod != 0, axis=0)
    # if overlapping, two diplotypes must share defining genotypes at every position and there cannot be an empty column
    status_sharing_definitions = np.all(non_empty_cols)

    # find diplotypes that is not equivalent to, inclusive of, or included in g1
    status_overlapping = set(g_subs.flatten()) == {1, 0, -1}

    # check whether each diplotype have exclusive defining genotypes
    status = status_sharing_definitions * status_overlapping

    return status


def find_pairwise_relationship(d1: str, d2: str, g1: np.ndarray, g2: np.ndarray) -> str:
    """
    a wrapper function to determine the relationship between two diplotypes
    the input d1 and d2 must share definitions

    :param d1: the name of diplotype 1
    :param d2: the name of diplotype 2
    :param g1: the numpy aray of diplotype 1 (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param g2: the numpy aray of diplotype 1 (M, P)
    :return: a string specifying the relationship between d1 and d2,
            including "equivalent", "overlapping", "included", and "inclusive"
    """
    # compare g1 with g2
    status_equivalent = is_equivalent(g1, g2)
    status_included = is_included(g1, g2)
    status_inclusive = is_inclusive(g1, g2)
    status_overlapping = is_overlapping(g1, g2)

    # g1 and g2 must have only one True status
    status_sum = sum([status_equivalent, status_included, status_inclusive, status_overlapping])
    if status_sum != 1:
        print(f'Something is wrong between {d1} and {d2}.')
        sys.exit(1)
    elif status_equivalent:
        relationship = 'equivalent'
    elif status_included:
        relationship = d1 + ' is included in ' + d2
    elif status_inclusive:
        relationship = d1 + ' is inclusive of ' + d2
    else:
        relationship = 'overlapping'

    return relationship


def find_possible_calls(g1: np.ndarray,
                        definition_arrays: np.ndarray,
                        diplotype_names: list[str]) -> list[str]:
    """
    find all alternative calls for one diplotype
    :param g1: the numpy aray of diplotype 1 (M, P)
                M = possible genotypes at a position for a diplotype
                P = allele-defining position for a gene
    :param definition_arrays: a numpy array of all diplotype definitions (D, M, P)
                D = the total number of diplotypes for a gene
    :param diplotype_names: a list of diplotype names
    :return: a list of diplotypes that share definitions with g1 and can show up as alternative calls for g1
    """
    # check whether g1 shares definitions with any other diplotype definitions
    status_sharing_definitions: np.ndarray = is_sharing_definitions(g1, definition_arrays)
    # check whether each of the diplotypes differs from g1
    status_discrepant: np.ndarray = is_discrepant(g1, definition_arrays)
    # sanity check to make sure no diplotypes can be discrepant but share definitions at the same time
    if (status_sharing_definitions == status_discrepant).any():
        print(f'Something went wrong. Some diplotypes are different but also the same at the same time.')
        sys.exit(1)
    # continue to the next diplotype if the current one is discrepant from all other diplotypes
    if all(status_discrepant):
        print(f'Something went wrong. A diplotype should at least share the allele definition with itself.')
        sys.exit(1)

    # get the names of the diplotypes that share definitions with g1
    diplotypes_sharing_definitions = [d for d, s in zip(diplotype_names, status_sharing_definitions) if s]

    return diplotypes_sharing_definitions


def find_all_possible_calls(diplotype_definitions: dict[str, dict[str, Set[str]]],
                            allele_defining_variants: dict[str, dict[str, str]]) -> dict:
    """
    find the alternative calls for all diplotypes.
    check out the function 'compare_diplotype_pairs' for how the relationships between diplotypes were determined.

    :param diplotype_definitions: a dictionary of diplotype definitions
    :param allele_defining_variants: a dictionary of all allele-defining variants for a gene
    :return: a dictionary of all the alternative calls one will find for a gene
        {diplotype: [list of all possible calls]} or
        {wobble diplotypes: {case123: [list of all possible calls]}} ## not all wobble diplotypes have multiple cases
    """
    # initialize an empty dictionary to store the alternative calls
    dict_possible_calls: dict = dict()

    # convert the diplotype definition dictionaries to numpy arrays of genotypes for each diplotype
    # vectorized computation by numpy arrays speeds up diplotype comparison
    definition_arrays: np.ndarray = convert_dictionary_to_numpy_array(diplotype_definitions,
                                                                      allele_defining_variants)

    # get the list of wobble diplotypes
    wobble_diplotypes: list[str] = []
    for d, v in diplotype_definitions.items():
        for g in v.values():
            if len(g) > 1:
                wobble_diplotypes.append(d)
                break

    # get the list of diplotypes that have overlapping definitions with other diplotype
    diplotype_names: list[str] = [*diplotype_definitions]
    # get the number of diplotypes
    n_diplotypes = len(diplotype_definitions)

    # loop over all diplotypes
    for i in range(n_diplotypes):

        # get the genotype array of the target diplotype
        d1 = diplotype_names[i]
        g1 = definition_arrays[i, :, :]

        # find the alternative calls of a diplotype
        possible_calls = find_possible_calls(g1, definition_arrays, diplotype_names)

        # if there is no other alternative calls except d1 itself, move on to the next diplotype
        if [d1] == possible_calls:
            continue
        # if d1 is not a wobble diplotype, add the alternative calls to the summary dictionary and continue
        elif d1 not in wobble_diplotypes:
            dict_possible_calls[d1] = possible_calls
        else:
            # add the list of alternative calls as a frozenset to the summary dictionary
            # set is not hashable, so use a frozenset instead
            # set is desired here so that the summary dictionary will not have redundant sets of alternative calls
            possible_call_sets: set[frozenset] = {frozenset(possible_calls)}

            for d2 in possible_calls:
                d2_idx = diplotype_names.index(d2)
                g2 = definition_arrays[d2_idx, :, :]

                # find the relationship between d1 and d2
                relationship = find_pairwise_relationship(d1, d2, g1, g2)

                # if d1 overlaps with d2, or d1 includes d2, find the more limited set of possible calls
                if relationship == 'overlapping' or 'inclusive' in relationship:
                    # find the intersecting definition of g1 and g2
                    g3 = g1 * g2

                    # find the alternative calls of g3
                    possible_calls = find_possible_calls(g3, definition_arrays, diplotype_names)

                    # add the alternative calls for g3 to the summary dictionary
                    # this is one of the results one will observe in the autogenerated tests
                    possible_call_sets.add(frozenset(possible_calls))

            # add the frozensets to the summary dictionary
            if len(possible_call_sets) == 1:
                dict_possible_calls[d1] = list(list(possible_call_sets)[0])
            else:
                # initialize a dictionary item for d1
                dict_possible_calls[d1] = dict()
                # add different sets of alternative calls to the summary dictionary for d1
                for j, x in enumerate(possible_call_sets):
                    dict_possible_calls[d1]['case_' + str(j+1)] = list(x)

    # return the summary dictionary of alternative calls for all diplotypes
    return dict_possible_calls


def get_actual_and_alternative_calls(possible_calls: list[str], diplotype_scores: dict[str, int]) \
                                     -> Tuple[str, str]:

    # initialize output variables
    actual_calls: list[str] = []
    alternative_calls: list[str] = []

    # get the scores for the possible diplotypes
    scores = np.array([diplotype_scores[d] for d in possible_calls])
    # get the diplotypes with the max scores
    max_idx = np.where(scores == np.max(scores))[0]

    # separate the actual and alternative calls based on diplotype scores
    for i, d in enumerate(possible_calls):
        (alternative_calls, actual_calls)[i in max_idx].append(d)

    # return str
    return ','.join(actual_calls), ','.join(alternative_calls)


def find_predicted_calls(dict_possible_calls: dict[str, str],
                         diplotype_scores: dict[str, int]) -> dict[list[str], list[str], list[str]]:
    # initialize the return dictionary
    predict_calls_cols = ['expected', 'actual', 'alternative']
    predicted_calls: dict = {key: [] for key in predict_calls_cols}

    for expected_call, possible_calls in dict_possible_calls.items():

        # if there is no multiple cases of possible calls
        if isinstance(possible_calls, list):
            # get the lists of actual calls and alternative calls based on scores
            actual_calls, alternative_calls = get_actual_and_alternative_calls(possible_calls, diplotype_scores)

            # add to the summary dictionary
            predicted_calls['expected'].append(expected_call)
            predicted_calls['actual'].append(actual_calls)
            predicted_calls['alternative'].append(alternative_calls)
        # some wobble diplotypes will have multiple cases of possible calls given different definitions
        elif isinstance(possible_calls, dict):
            # iterate over different sets of possible calls for a wobble diplotype
            for case in possible_calls.values():

                # get the lists of actual calls and alternative calls based on scores
                actual_calls, alternative_calls = get_actual_and_alternative_calls(case, diplotype_scores)

                # add to the summary dictionary
                predicted_calls['expected'].append(expected_call)
                predicted_calls['actual'].append(actual_calls)
                predicted_calls['alternative'].append(alternative_calls)
        else:
            print(f'Something went wrong for {expected_call}. Data type was wrong for possible calls.')

    # return the expected, actual, and alternative calls
    return predicted_calls


def compare_diplotype_pairs(diplotype_definitions: dict,
                            allele_defining_variants: dict[str, dict[str, str]]) -> dict:
    """
    compare unphased diplotype definitions
    G = a MxP matrix where rows are all possible nucleotide pairs in a gene (M) and columns are positions (P)

    Possible outcomes of diplotype comparison
    1. Sharing definitions
        - between any diplotypes (G1, G2)
        - criteria:
            1) G1 * G2 = G3. G3 cannot have empty columns.
                empty column = G1 and G2 do not share any definitions at a position

        Four possible relationships between diplotypes that share definitions
        1. equivalent
            - between any diplotypes (G1, G2)
            - criteria:
                1) G1 shares a definition with G2
                2) G1 - G2 = G3 = {0}
            - alternative criteria:
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 = {0}

        2. included
            - G1 (a wobble or a non-wobble diplotype) is included in G2 (another wobble diplotype) as if a sub-allele
            - criteria:
                1) G1 shares a definition with G2
                2) G1 != G2 (not equivalent)
                3) G1 <= G2
            - alternative criteria
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 where G4 must be made up by {0, -1}

            criteria:
            1) (integer operation) G1 != G2 and G1 <= G2
            1) (alternative) G1 - G2 = G3 with only {0,-1}. G3 cannot have empty columns.

        3. inclusive
            - G1 (a wobble diplotype) includes G2 (other wobble and non-wobble diplotypes)
            - criteria:
                1) G1 shares a definition with G2
                2) G1 != G2 (not equivalent)
                3) G1 >= G2
            - alternative criteria
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 where G4 must be made up by {0, -1}

        4. overlapping
            - G1 (a wobble diplotype) overlaps with G2 (another wobble diplotype)
            - criteria:
                1) G1 shares a definition with G2
                2) G1 != G2 (not equivalent)
                3) G1 is not greater than, nor less than, G2
            - alternative criteria
                1) G1 * G2 = G3. G3 cannot have empty columns.
                2) G1 - G2 = G4 where G4 must be made up by {1, 0, -1}

    2. discrepant
        - between any diplotypes (G1, G2)
        - criteria:
            1) G1 - G2 = G3 where G3 must be made up by {1, 0, -1}
            2) G1 * G2 = G4. G4 has empty columns
                empty column = G1 and G2 do not share any defining genotypes at a position


    :param allele_defining_variants:
    :param diplotype_definitions:
    :return: a dictionary {
        'diplotype_1': {
            score: 123,
            diplotype_name: relation
        } }
    """
    # set up an empty dictionary to store diplotype comparison results
    pairwise_comparison = dict()

    # convert the diplotype definition dictionaries to numpy arrays of genotypes for each diplotype
    # vectorized computation by numpy arrays speeds up diplotype comparison
    definition_arrays: np.ndarray = convert_dictionary_to_numpy_array(diplotype_definitions,
                                                                      allele_defining_variants)

    # get the list of diplotype names
    diplotype_names: list[str] = [*diplotype_definitions]
    # get the number of diplotypes
    n_diplotypes = len(diplotype_definitions)

    # iterate over diplotypes
    for i in range(n_diplotypes):

        # get the genotype array of the target diplotype
        d1 = diplotype_names[i]
        g1 = definition_arrays[i, :, :]

        # check whether g1 shares definitions with any other diplotype definitions
        status_sharing_definitions = is_sharing_definitions(g1, definition_arrays)
        # check whether each of the diplotypes differs from g1
        status_discrepant = is_discrepant(g1, definition_arrays)
        # sanity check to make sure no diplotypes can be discrepant but share definitions at the same time
        if any([status_sharing_definitions == status_discrepant]):
            print(f'Something went wrong. Some diplotypes are different but also the same at the same time.')
            sys.exit(1)
        # continue to the next diplotype if the current one is discrepant from all other diplotypes
        if all(status_discrepant):
            print(f'Something went wrong. A diplotype should at least share the allele definition with itself.')
            sys.exit(1)

        # specify the relationship between diplotypes that share definitions
        for j, s in enumerate(status_sharing_definitions):
            if i == j:  # skip self
                continue
            elif s:  # s == True meaning two diplotypes share definitions
                # get the definition array of another diplotype
                d2 = diplotype_names[j]
                g2 = definition_arrays[j, :, :]

                # find the relationship between d1 and d2
                relationship = find_pairwise_relationship(d1, d2, g1, g2)

                # add diplotype and its score to the output dictionary
                if d1 not in pairwise_comparison:
                    pairwise_comparison[d1] = {}

                # add the relationship between g1 and g2 to a dictionary
                pairwise_comparison[d1][d2] = relationship

            else:
                continue

    return pairwise_comparison
