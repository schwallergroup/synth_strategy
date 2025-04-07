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
    This function detects the use of TMS-protected fragments in the synthesis.
    """
    has_tms_protected_fragment = False

    def dfs_traverse(node):
        nonlocal has_tms_protected_fragment

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for TMS group in reactants
                tms_pattern = Chem.MolFromSmarts("[#14]([#6])([#6])[#6]")

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(tms_pattern):
                            print("Detected TMS-protected fragment")
                            has_tms_protected_fragment = True
                            break
                    except:
                        continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_tms_protected_fragment
