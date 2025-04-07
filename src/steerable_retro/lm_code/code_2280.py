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
    This function detects sequential functionalization of a triazine core
    with at least two substitution steps.
    """
    triazine_substitution_count = 0
    triazine_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#7][#6][#7]1")

    def dfs_traverse(node):
        nonlocal triazine_substitution_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a nucleophilic aromatic substitution on triazine
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(triazine_pattern):
                    # Look for chlorine replacement on triazine
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(triazine_pattern):
                            # Check if reactant has Cl attached to triazine
                            cl_triazine_pattern = Chem.MolFromSmarts(
                                "[Cl]-[#6]1[#7][#6][#7][#6][#7]1"
                            )
                            if reactant_mol.HasSubstructMatch(cl_triazine_pattern):
                                triazine_substitution_count += 1
                                print(f"Found triazine substitution: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found at least 2 sequential triazine substitutions
    result = triazine_substitution_count >= 2
    print(
        f"Triazine sequential functionalization detected: {result} (count: {triazine_substitution_count})"
    )
    return result
