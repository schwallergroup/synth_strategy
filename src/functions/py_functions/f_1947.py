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
    This function detects the presence of a piperazine scaffold that undergoes
    multiple substitution reactions at different stages of the synthesis.
    """
    piperazine_pattern = Chem.MolFromSmarts("[N]1[C][C][N][C][C]1")
    substitution_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check if product contains piperazine
                if product_mol and product_mol.HasSubstructMatch(piperazine_pattern):
                    # Check if this is a substitution reaction on piperazine
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            piperazine_pattern
                        ):
                            # This is potentially a substitution on piperazine
                            substitution_depths.append(depth)
                            print(f"Piperazine substitution detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if piperazine undergoes multiple substitutions
    return len(substitution_depths) >= 2
