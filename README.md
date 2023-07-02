# Tanimoto Similarity

## Introduction

The tanimoto.py cli program calculates the Tanimoto similarity between two lists of SMILES strings. The Tanimoto similarity is a measure of the similarity between two sets. In this case, the two sets are the sets of substructures in each SMILES string. The Tanimoto similarity is calculated as follows:

$$
T(A,B) = \frac{|A \cap B|}{|A \cup B|}
$$

where $A$ and $B$ are the two sets, and $|A|$ is the number of elements in $A$. The Tanimoto similarity is a number between 0 and 1, where 0 means that the two sets are completely different, and 1 means that the two sets are identical.


## Installation

The tanimoto.py cli program requires Python 3.6 or higher. It also requires the RDKit cheminformatics software. The RDKit can be installed using conda:

```
conda install -c rdkit rdkit
```

The tanimoto.py cli program can be installed by cloning the GitHub repository:

```
git clone
```

## Usage

The tanimoto.py cli program can be run from the command line as follows:

```
python tanimoto.py --file1 <filename> --file2 <filename>
```

The program takes two arguments, the names of two files containing SMILES strings. The program calculates the Tanimoto similarity between each pair of SMILES strings, and writes the results to a file called `results.csv`.

If you want to specify the number of most similar molecules to return, you can use the `--top` argument:

```
python tanimoto.py --file1 <filename> --file2 <filename> --top 10
```

## Credits

This project was developed by <name> as part of the <name> course at <institution>.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
