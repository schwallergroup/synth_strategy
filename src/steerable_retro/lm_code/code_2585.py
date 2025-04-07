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
    Detects if the synthesis uses a thiourea intermediate to form a heterocycle.
    This is a more general version of the thiazole formation strategy.
    """
    # Track thiourea presence and subsequent use
    thiourea_found = False
    thiourea_used_for_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal thiourea_found, thiourea_used_for_heterocycle

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for thiourea in reactants
            thiourea_pattern = Chem.MolFromSmarts("[#6](=[#16])[#7]")
            thiourea_in_reactants = False

            for reactant_smiles in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                if reactant_mol and reactant_mol.HasSubstructMatch(thiourea_pattern):
                    thiourea_in_reactants = True
                    thiourea_found = True
                    break

            # Check for heterocycle formation from thiourea
            if thiourea_in_reactants:
                product_mol = Chem.MolFromSmiles(product_smiles)
                # Check for common heterocycles
                thiazole_pattern = Chem.MolFromSmarts("[s]1[c][n][c][c]1")
                oxazole_pattern = Chem.MolFromSmarts("[o]1[c][n][c][c]1")
                imidazole_pattern = Chem.MolFromSmarts("[n]1[c][n][c][c]1")

                if product_mol and (
                    product_mol.HasSubstructMatch(thiazole_pattern)
                    or product_mol.HasSubstructMatch(oxazole_pattern)
                    or product_mol.HasSubstructMatch(imidazole_pattern)
                ):
                    thiourea_used_for_heterocycle = True
                    print(f"Found heterocycle formation from thiourea at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if thiourea_found and thiourea_used_for_heterocycle:
        print("Detected thiourea-mediated heterocycle synthesis")
        return True
    return False
