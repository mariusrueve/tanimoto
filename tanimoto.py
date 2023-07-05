"""
This script calculates the Tanimoto similarity between two lists of molecules.
The input files should contain a list of SMILES strings, one per line.
The output is a csv file with the top n most similar molecules from mol_list2 for each molecule in mol_list1.
If --top is not specified, the top 10 most similar molecules are returned.
"""

import argparse
import os
import pathlib as pl
import sys
import time

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
    parser.add_argument(
        "--top",
        type=int,
        help="Number of most similar molecules to return",
        required=False,
    )
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
    """Returns a list with mol1, mol2, similarity"""
    similarity_list = []
    for mol1 in mol_list1:
        for mol2 in mol_list2:
            if mol1 is None or mol2 is None:
                continue
            else:
                similarity = get_similarity(mol1, mol2)
            similarity_list.append(
                [Chem.MolToSmiles(mol1), Chem.MolToSmiles(mol2), similarity]
            )
    return similarity_list


def get_top_n(similarity_list, n):
    """Returns the top n most similar molecules from mol_list2 for each molecule in mol_list1"""
    df = pd.DataFrame(similarity_list, columns=["mol1", "mol2", "similarity"])
    df["mol1"] = pd.Categorical(df["mol1"], categories=df["mol1"].unique())
    df = df.sort_values(["mol1", "similarity"], ascending=[True, False])
    df = df.groupby("mol1").head(n)
    return df


def main():
    args = get_arguments()
    mol_list1 = get_mol_list(
        os.path.join(pl.Path(__file__).parent.absolute(), args.file1)
    )
    mol_list2 = get_mol_list(
        os.path.join(pl.Path(__file__).parent.absolute(), args.file2)
    )
    similarity_list = get_similarity_list(mol_list1, mol_list2)
    # call get_top_n with the number of most similar molecules to return.
    # if --top is not specified, return the top 10 most similar molecules
    if args.top:
        df = get_top_n(similarity_list, args.top)
    else:
        df = get_top_n(similarity_list, 10)

    output_location = os.path.join(pl.Path(__file__).parent.absolute(), "results.csv")
    if args.output:
        output_location = os.path.join(pl.Path(__file__).parent.absolute(), args.output)

    df.to_csv(output_location, index=False)
    print(f"Results written to '{output_location}'", file=sys.stderr)


if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(f"Time elapsed: {end - start} seconds", file=sys.stderr)
