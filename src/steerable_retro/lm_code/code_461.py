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
    Detects if the synthetic route involves multiple C-N bond formations (at least 3).
    """
    c_n_bond_formations = 0

    def dfs_traverse(node):
        nonlocal c_n_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Create molecules from SMILES
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check for C-N bond formation
                    c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                    reactants_c_n_count = len(reactants_mol.GetSubstructMatches(c_n_pattern))
                    product_c_n_count = len(product_mol.GetSubstructMatches(c_n_pattern))

                    if product_c_n_count > reactants_c_n_count:
                        c_n_bond_formations += 1
                        print(f"Found C-N bond formation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return c_n_bond_formations >= 3
