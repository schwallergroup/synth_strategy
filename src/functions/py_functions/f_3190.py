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
    This function detects if the synthetic route involves indole ring formation.
    """
    indole_formed = False
    indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")

    def dfs_traverse(node):
        nonlocal indole_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains indole but reactants don't
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(indole_pattern):
                    reactants_have_indole = False
                    for r_smi in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if r_mol and r_mol.HasSubstructMatch(indole_pattern):
                            reactants_have_indole = True
                            break

                    if not reactants_have_indole:
                        print(f"Indole formation detected in reaction: {rsmi}")
                        indole_formed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return indole_formed
