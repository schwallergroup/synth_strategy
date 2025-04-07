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
    This function detects if the synthesis involves N-alkylation of a pyrazole ring.
    """
    has_pyrazole_alkylation = False

    def dfs_traverse(node):
        nonlocal has_pyrazole_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(reactants_mols):
                # Check for pyrazole in reactants
                pyrazole_pattern = Chem.MolFromSmarts("[n]1[n]cc[c]1")

                # Check for N-alkylation: pyrazole-NH to pyrazole-N-C
                pyrazole_nh_pattern = Chem.MolFromSmarts("[nH]1[n]cc[c]1")
                pyrazole_nc_pattern = Chem.MolFromSmarts("[n]1([C])[n]cc[c]1")

                if any(
                    mol.HasSubstructMatch(pyrazole_nh_pattern) for mol in reactants_mols
                ) and product_mol.HasSubstructMatch(pyrazole_nc_pattern):
                    has_pyrazole_alkylation = True
                    print(f"Detected pyrazole N-alkylation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_pyrazole_alkylation
