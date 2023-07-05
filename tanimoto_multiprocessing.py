"""Write a program that takes two lists of molecules as input and calculate the tanimoto similarity using the morgan fingerprint between each molecule in file 1
and each molecule in file 2. Use multiprocessing to speed up the calculation."""

import argparse
import os
import pathlib as pl
import sys
import time

import multiprocessing as mp
from multiprocessing import Pool

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Tanimoto similarity between two lists of molecules"
    )
    parser.add_argument("--file1", type=str, help="Path to file1", required=True)
    parser.add_argument("--file2", type=str, help="Path to file2", required=True)
    # optional argument to specify location of output file
    parser.add_argument(
        "--output", type=str, help="Path to output file", required=False
    )

    args = parser.parse_args()
    return args


def get_mol_list(file: os.path):
    """Reads a file with SMILES strings and returns a list of RDKit molecules"""
    with open(file) as f:
        smiles = f.readlines()
    smiles = [x.strip() for x in smiles]
    mol_list = [Chem.MolFromSmiles(x) for x in smiles]
    # get relative path of file
    print(
        f"Succesfully read '{file}', Number of molecules: {len(mol_list)}",
        file=sys.stderr,
    )
    return mol_list


def get_similarity(mol1, mol2):
    """Calculates the Tanimoto similarity between two molecules"""
    fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
    fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
    # check if both fingerprints are valid
    if fp1 is None or fp2 is None:
        similarity = 0
    else:
        similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity


def get_similarity_list(mol_list1, mol_list2):
    """Calculates the Tanimoto similarity between each molecule in mol_list1 and each molecule in mol_list2 while using multiprocessing"""
    # create a list of tuples with all combinations of molecules
    mol_combinations = [(x, y) for x in mol_list1 for y in mol_list2]
    # wirte the line above not as a one liner
    mol_combinations = []
    for x in mol_list1:
        for y in mol_list2:
            if x is None or y is None:
                continue
            else:
                mol_combinations.append((x, y))

    # create a pool of workers
    pool = Pool(mp.cpu_count())
    # calculate the similarity for each tuple in mol_combinations
    print(mol_combinations)
    similarity_list = pool.starmap(get_similarity, mol_combinations)
    pool.close()
    pool.join()
    return similarity_list


def get_top_matches(similarity_list, mol_list1, mol_list2):
    """Returns a dataframe with the best single match for each molecule in mol_list1"""
    # create a list of tuples with all combinations of molecules
    mol_combinations = []
    for x in mol_list1:
        for y in mol_list2:
            if x is None or y is None:
                continue
            else:
                mol_combinations.append((Chem.MolToSmiles(x), Chem.MolToSmiles(y)))
    # create a dataframe with the results
    df = pd.DataFrame(mol_combinations, columns=["mol1", "mol2"])
    df["similarity"] = similarity_list
    # Manipulate the dataframe such that it only contains the best match for each molecule in mol_list1
    df = df.sort_values(by=["similarity"], ascending=False)
    df_top = df.drop_duplicates(subset=["mol1"], keep="first")
    return df_top


def main():
    args = get_arguments()
    file1 = args.file1
    file2 = args.file2
    output = args.output

    # get the molecules from the files
    mol_list1 = get_mol_list(file1)
    mol_list2 = get_mol_list(file2)

    # calculate the similarity between all molecules
    similarity_list = get_similarity_list(mol_list1, mol_list2)

    # get the top matches
    df_top = get_top_matches(similarity_list, mol_list1, mol_list2)

    # write the results to a file
    file_name = pl.Path(file1).stem + "_" + pl.Path(file2).stem + ".csv"
    if output:
        file_name = output
    df_top.to_csv(file_name, index=False)
    print(f"Results written to '{file_name}'", file=sys.stderr)


if __name__ == "__main__":
    time_start = time.time()
    main()
    time_end = time.time()
    print(f"Runtime: {time_end - time_start:.2f} seconds", file=sys.stderr)
