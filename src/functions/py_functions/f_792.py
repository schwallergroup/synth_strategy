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
    Detects if the synthesis route involves protection of a phenol as a methoxy group.
    """
    protection_found = False

    def dfs_traverse(node):
        nonlocal protection_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for phenol to methoxy transformation
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and any(r for r in reactant_mols if r):
                # Check for phenol pattern in reactants
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                # Check for methoxy pattern in product
                methoxy_pattern = Chem.MolFromSmarts("[c][O][CH3]")

                has_phenol = any(
                    r.HasSubstructMatch(phenol_pattern) for r in reactant_mols if r
                )
                has_methoxy = product_mol.HasSubstructMatch(methoxy_pattern)

                if has_phenol and has_methoxy:
                    protection_found = True
                    print("Phenol protection found")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Phenol protection strategy: {protection_found}")
    return protection_found
