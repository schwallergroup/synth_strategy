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
    This function detects the use of a Boc-protected amine as a nucleophile in the synthesis.
    """
    boc_protected_amine = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protected_amine

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for Boc-protected amine in reactants
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        boc_pattern = Chem.MolFromSmarts("[CH3]C([CH3])([CH3])[O]C(=O)[NH]")
                        if reactant_mol.HasSubstructMatch(boc_pattern):
                            boc_protected_amine = True
                            print(f"Detected Boc-protected amine at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return boc_protected_amine
