import argparse
import os
import sys
import time

import multiprocessing as mp
from multiprocessing import Pool

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# q: How can I make this script more efficient?
# a: Use a fingerprint list instead of a molecule list

# q: Doesn't i already use a fingerprint list?
# a: Yes, but it is calculated on the fly. Instead, we can calculate the fingerprints for each molecule in the list once and store them in a list.

# q: How do I do that?
# a: See the calculate_fingerprint function below. It takes a list of molecules as input and returns a list of fingerprints.
# The fingerprints are calculated using the RDKit function GetMorganFingerprintAsBitVect. The function takes three arguments:
# the molecule, the radius of the fingerprint and the number of bits.
# The radius is the number of bonds that are considered when calculating the fingerprint.
# The number of bits is the length of the fingerprint. The default value for the radius is 2 and the default value for the number of bits is 2048.
# In this case, we use a radius of 2 and a number of bits of 1024. The radius and number of bits are chosen to be the same as in the original script.


def get_arguments():
    # parse argument for file1
    parser = argparse.ArgumentParser(
        description="Tanimoto similarity between two lists of molecules"
    )
    parser.add_argument(
        "--file1",
        type=str,
        help="Path to file1",
        required=True,
    )
    parser.add_argument(
        "--file2",
        type=str,
        help="Path to file2",
        required=True,
    )

    args = parser.parse_args()
    return args


def read_mol_list(file: os.path):
    """Reads a file with SMILES strings and returns a list of RDKit molecules"""
    with open(file) as f:
        smiles = f.readlines()
    smiles = [x.strip() for x in smiles]
    mol_list = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            mol_list.append(mol)
    # get relative path of file
    print(
        f"Succesfully read '{file}', Number of molecules: {len(mol_list)}",
        file=sys.stderr,
    )
    return mol_list


def calculate_fingerprint(mol_list):
    """Calculates the Morgan fingerprint for each molecule in a list"""
    fingerprint_list = [
        rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        for mol in mol_list
    ]

    return fingerprint_list


def get_similarity(fingerprint1, fingerprint2):
    """Calculates the Tanimoto similarity between two molecules"""
    similarity = DataStructs.TanimotoSimilarity(fingerprint1, fingerprint2)
    return similarity


def main():
    args = get_arguments()
    file1 = args.file1
    file2 = args.file2

    # read molecules from file
    mol_list1 = read_mol_list(file1)
    mol_list2 = read_mol_list(file2)

    # calculate fingerprints
    fingerprint_list1 = calculate_fingerprint(mol_list1)
    fingerprint_list2 = calculate_fingerprint(mol_list2)

    similarity_list = []
    for i, fingpergrint1 in enumerate(fingerprint_list1):
        best_similarity = 0
        for j, fingerprint2 in enumerate(fingerprint_list2):
            similarity = get_similarity(fingpergrint1, fingerprint2)
            if similarity > best_similarity:
                best_similarity = similarity
                best_mol_index = j
        similarity_list.append((i, best_mol_index, best_similarity))

    # sort list by similarity
    similarity_list.sort(key=lambda x: x[2], reverse=True)

    # write each similarity to file
    for i, j, similarity in similarity_list:
        print(
            f"{Chem.MolToSmiles(mol_list1[i])}\t{Chem.MolToSmiles(mol_list2[j])}\t{similarity}"
        )


if __name__ == "__main__":
    start_time = time.time()
    main()
    print(f"--- {time.time() - start_time} seconds ---", file=sys.stderr)
