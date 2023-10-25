#!/usr/bin/env python3

__author__ = 'BinglanLi'

from glob import glob
import sys
from pathlib import Path
import pandas as pd

import utilities as util

if __name__ == "__main__":
    import argparse

    # describe the tool
    parser = argparse.ArgumentParser(description='Identify overlapping diplotypes in unphased genetic data.')

    # input arguments
    parser.add_argument("-i", "--input",
                        type=str,
                        required=False,
                        # default='/Users/binglan/Documents/GitHub/PharmCAT/src/main/resources/org/pharmgkb/pharmcat'
                        #         '/definition/alleles/',
                        help="File name or the directory to the PharmCAT JSONs of allele definitions.")
    # input arguments
    parser.add_argument("-m", "--missing-pos",
                        type=str,
                        required=False,
                        metavar='<file>',
                        default='',
                        help="List of missing positions. Each row is a set of missing positions.")

    # output args
    parser.add_argument("-o", "--output-dir",
                        type=str,
                        required=False,
                        metavar='<dir>',
                        help="(Optional) directory for outputs.  Defaults to the current working directory.")
    parser.add_argument("-bf", "--base-filename",
                        type=str,
                        required=False,
                        default='predicted_pharmcat_calls',
                        help="(Optional) base filename of the output file.")

    # parse arguments
    args = parser.parse_args()

    try:
        # define output file name
        m_output_basename: str = args.base_filename
        # define output file path
        m_output_dir: Path = args.output_dir if args.output_dir else Path.cwd()
        if not args.output_dir:
            print(f'Output directory is not specified.\n'
                  f'\tUse the current working directory: {m_output_dir}')

        # define input
        m_input: Path = Path(args.input).expanduser() if args.input else Path.cwd()
        allele_definition_jsons: list = []
        if not args.input:
            print(f'No input was provided. '
                  f'Looking for the allele definition JSON files under the current working directory: {m_input}')
        elif m_input.is_file():
            allele_definition_jsons: list = [m_input]
        elif m_input.is_dir():
            print(f'Looking for the allele definition JSON files under {m_input}')
            # get the list allele definition files
            allele_definition_jsons: list = glob(str(m_input.joinpath('*_translation.json')))
        # check whether there is any json files
        if len(allele_definition_jsons) == 0:
            print(f'No allele definition JSON files are found under {m_input}')
            sys.exit(1)

        # process each gene
        for json_file in allele_definition_jsons:
            # read the json file
            print(f'Processing {json_file}')
            json_data: dict = util.read_json(json_file)

            # get necessary information of the allele-defining positions, including hgvs names, positions, references
            allele_defining_variants: dict = util.get_allele_defining_variants(json_data)

            # read in the list of alleles, allele-defining positions, and defining genotypes
            allele_definitions = util.get_allele_definitions(json_data, allele_defining_variants)

            # update the allele_definitions and fill empty cells with reference genotypes
            allele_definitions_complete: dict = util.fill_definitions_with_references(allele_definitions)

            # get a dictionary of diplotype definitions for comparison
            # {diplotype_name: {position: {unique combinations of genotypes at a position} } }
            diplotype_definitions: dict = util.get_diplotype_definition_dictionary(
                allele_definitions_complete, allele_defining_variants)

            # find all possible outcomes of alternative calls for each diplotype
            dict_possible_calls: dict = util.find_all_possible_calls(diplotype_definitions,
                                                                     allele_defining_variants)

            # get a dictionary of calculated diplotype scores based on the number of core allele-defining positions
            # first, calculate the haplotype scores
            allele_scores: dict[str, int] = util.count_allele_scores(allele_definitions)
            # then, calculate the diplotype scores by adding up the allele scores
            diplotype_scores: dict[str, int] = util.count_diplotype_scores(allele_scores)
            # summarize expected vs actual calls
            dict_predicted_calls: dict = util.find_predicted_calls(dict_possible_calls, diplotype_scores)

            # convert the python dictionary to a pandas data frame for output
            df_predicted_calls = pd.DataFrame.from_dict(dict_predicted_calls)
            # get the gene name
            gene: str = json_file.split('/')[-1].split('_')[0]
            # add gene name to the data frame for readability in the output
            df_predicted_calls['gene'] = gene
            # specify output name for a single gene
            m_output_basename_gene = m_output_basename + "_" + gene + ".tsv"
            # write to an output
            df_predicted_calls.to_csv(m_output_basename_gene, sep="\t", index=False)

            # compare diplotype pairs and identify pairwise relationships
            # pairwise_comparisons: dict = util.compare_diplotype_pairs(diplotype_definitions, allele_defining_variants)

    except Exception as e:
        print(e)
        sys.exit(1)
