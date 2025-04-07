#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects a strategy involving functionalization of an imidazole core
    at multiple positions.
    """
    imidazole_functionalizations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for imidazole core
            imidazole_pattern = Chem.MolFromSmarts("c1ncnc1")

            # Check if product has imidazole
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(imidazole_pattern):
                # This reaction involves an imidazole

                # Check for different types of functionalizations
                n_alkylation_pattern = Chem.MolFromSmarts("[n][CH2]c")
                c_oxidation_pattern = Chem.MolFromSmarts("c[CH]=O")

                if product_mol.HasSubstructMatch(n_alkylation_pattern):
                    imidazole_functionalizations.append(("n_alkylation", depth))
                    print(f"Detected imidazole N-alkylation at depth {depth}")

                if product_mol.HasSubstructMatch(c_oxidation_pattern):
                    imidazole_functionalizations.append(("c_oxidation", depth))
                    print(f"Detected imidazole C-oxidation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have at least 2 different types of functionalizations
    functionalization_types = set(f[0] for f in imidazole_functionalizations)
    multiple_functionalizations = len(functionalization_types) >= 2

    if multiple_functionalizations:
        print(f"Detected multiple types of imidazole functionalizations: {functionalization_types}")

    return multiple_functionalizations
