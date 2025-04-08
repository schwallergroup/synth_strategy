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
    Detects a synthetic strategy involving diarylamine formation through SNAr.
    """
    diarylamine_formed = False

    def dfs_traverse(node):
        nonlocal diarylamine_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for diarylamine formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    diarylamine_pattern = Chem.MolFromSmarts("c-[NH]-c")
                    if product_mol.HasSubstructMatch(diarylamine_pattern):
                        # Verify it's a new formation by checking reactants
                        has_pattern_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(diarylamine_pattern):
                                has_pattern_in_reactants = True
                                break

                        if not has_pattern_in_reactants:
                            diarylamine_formed = True
                            print("Diarylamine formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Diarylamine formation strategy detected: {diarylamine_formed}")
    return diarylamine_formed
