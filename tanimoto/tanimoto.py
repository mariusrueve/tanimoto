import rdkit
import numpy as np
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

# suppress rdkit warnings
lg = rdkit.RDLogger.logger()
lg.setLevel(rdkit.RDLogger.CRITICAL)



def read_molecules(path):
    mols = Chem.SmilesMolSupplier(path)
    mols = [mol for mol in mols if mol is not None]
    return mols


def tanimoto_coefficient(query_fingerprint, database_fingerprints):
    """Looks for the best match in the database for the query molecule."""
    # Calculate the Tanimoto coefficient
    similarities = [
        DataStructs.FingerprintSimilarity(query_fingerprint, fingerprint)
        for fingerprint in database_fingerprints
    ]

    # Return the best match
    return np.max(similarities)


def tanimoto_coefficient_list(query_file, database_file):
    """Looks for the best match in the database for each query molecule."""
    query_molecules = read_molecules(query_file)
    database_molecules = read_molecules(database_file)

    fpgen = AllChem.GetRDKitFPGenerator()
    query_fps = [fpgen.GetFingerprint(mol) for mol in query_molecules]
    database_fps = [fpgen.GetFingerprint(mol) for mol in database_molecules]

    similarities = [
        tanimoto_coefficient(query_fp, database_fps)
        for query_fp in query_fps
    ]

    return query_molecules, database_molecules, similarities


def print_results(query_molecules, database_molecules, similarities):
    if len(query_molecules) == len(database_molecules) == len(similarities):
        for query_molecule, database_molecule, similarity in zip(query_molecules, database_molecules, similarities):
            print(f'{Chem.MolToSmiles(query_molecule)},{Chem.MolToSmiles(database_molecule)},{similarity}')
