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
    Detects a synthetic strategy involving diaryl ether formation.
    """
    diaryl_ether_formed = False

    def dfs_traverse(node):
        nonlocal diaryl_ether_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for diaryl ether formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    diaryl_ether_pattern = Chem.MolFromSmarts("c-O-[CH2]-c")
                    if product_mol.HasSubstructMatch(diaryl_ether_pattern):
                        # Verify it's a new formation by checking reactants
                        has_pattern_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                diaryl_ether_pattern
                            ):
                                has_pattern_in_reactants = True
                                break

                        if not has_pattern_in_reactants:
                            diaryl_ether_formed = True
                            print("Diaryl ether formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Diaryl ether formation strategy detected: {diaryl_ether_formed}")
    return diaryl_ether_formed
