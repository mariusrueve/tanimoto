import rdkit
import numpy as np
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

# suppress rdkit warnings
lg = rdkit.RDLogger.logger()
lg.setLevel(rdkit.RDLogger.CRITICAL)


def read_molecules(path, verbose=False):
    mols = Chem.SmilesMolSupplier(path)
    initial_length = len(mols)
    mols = [mol for mol in mols if mol is not None]
    if verbose:
        logging.info(f"Number of molecules that could not be read: {initial_length - len(mols)}")
    return mols


def tanimoto_coefficient(query_fingerprint, database_fingerprints):
    """Looks for the best match in the database for the query molecule."""
    # Calculate the Tanimoto coefficient
    similarities = [
        DataStructs.FingerprintSimilarity(query_fingerprint, fingerprint)
        for fingerprint in database_fingerprints
    ]

    # get the index of the best match
    best_match_index = np.argmax(similarities)

    return best_match_index, similarities[best_match_index]


def tanimoto_coefficient_list(query_file, database_file, verbose=False):
    """Looks for the best match in the database for each query molecule."""
    query_molecules = read_molecules(query_file, verbose=verbose)
    if verbose:
        logging.info(f"Number of query molecules read: {len(query_molecules)}")
    database_molecules = read_molecules(database_file, verbose=verbose)
    if verbose:
        logging.info(f"Number of database molecules read: {len(database_molecules)}")

    fpgen = AllChem.GetRDKitFPGenerator()
    query_fps = [fpgen.GetFingerprint(mol) for mol in query_molecules]
    database_fps = [fpgen.GetFingerprint(mol) for mol in database_molecules]

    best_match_database_molecules = []
    similarities = []

    for query_fp in query_fps:
        best_match_index, similarity = tanimoto_coefficient(query_fp, database_fps)
        best_match_database_molecules.append(database_molecules[best_match_index])
        similarities.append(similarity)

    return query_molecules, best_match_database_molecules, similarities


def print_results(query_molecules, database_molecules, similarities):
    if len(query_molecules) == len(similarities):
        for query_molecule, database_molecule, similarity in zip(
            query_molecules, database_molecules, similarities
        ):
            print(
                f"{Chem.MolToSmiles(query_molecule)},{Chem.MolToSmiles(database_molecule)},{similarity}"
            )
