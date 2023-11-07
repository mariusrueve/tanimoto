
import argparse
import logging

from tanimoto.tanimoto import tanimoto_coefficient_list, print_results


def main():
    parser = argparse.ArgumentParser(description='Tanimoto coefficient calculator')
    parser.add_argument('query_molecules', type=str, help='path to a file that contains molecules that will be used as a query')
    parser.add_argument('database_molecules', type=str, help='path to a file that contains molecules that will be used as a database')
    parser.add_argument('-v','--verbose', action='store_true', help='enable logging')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    query_molecules, database_molecules, similarities = tanimoto_coefficient_list(args.query_molecules, args.database_molecules)
    
    print_results(query_molecules, database_molecules, similarities)


if __name__ == '__main__':
    main()