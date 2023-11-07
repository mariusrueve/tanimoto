# Tanimoto Similarity

## Introduction

The tanimoto.py cli program calculates the Tanimoto similarity between two lists of SMILES strings. The Tanimoto similarity is a measure of the similarity between two sets. In this case, the two sets are the sets of substructures in each SMILES string. The Tanimoto similarity is calculated as follows:

$$
T(A,B) = \frac{|A \cap B|}{|A \cup B|}
$$

where $A$ and $B$ are the two sets, and $|A|$ is the number of elements in $A$. The Tanimoto similarity is a number between 0 and 1, where 0 means that the two sets are completely different, and 1 means that the two sets are identical.

## Installation

The tanimoto.py program requires Python 3.6 or higher. The program can be installed using pip:

```bash
git clone https://github.com/mariusrueve/tanimoto.git
pip install tanimoto
```

## Usage

The tanimoto.py program can be used as a command line interface (CLI) program. The program takes two files as input, and calculates the Tanimoto similarity between the two files. The files should be in SMI format, with one SMILES string per line. The program can be run as follows:

```bash
tanimoto-cli file1.smi file2.smi
```

The program then prints the Tanimoto similarity between the two files to the terminal. If you want to save the output to a file, you can pipe the output to a file:

```bash
tanimoto-cli file1.smi file2.smi > output.csv
```


## Development

The tanimoto.py program is developed using Python 3.6. The program uses the following libraries:

- RDKit
- Num Py

