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
    Detects if the synthesis follows a linear pathway rather than convergent.
    """

    def is_linear_synthesis(node):
        # For reaction nodes, check if there are at most 2 reactants
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]
            reactants = rxn_smiles.split(">")[0].split(".")

            # Filter out common solvents and reagents that don't indicate convergent synthesis
            common_solvents = [
                "CC(=O)O",
                "C1CCOC1",
                "BrBr",
                "CC(C)O",
                "CCO",
                "CCOC",
                "O",
                "CO",
                "CC(C)(C)O",
            ]
            actual_reactants = [r for r in reactants if r not in common_solvents]

            if len(actual_reactants) > 2:
                print(f"Non-linear reaction (>2 reactants): {rxn_smiles}")
                return False

        # For molecule nodes, check if they have more than one reaction child
        if node["type"] == "mol" and not node.get("in_stock", False):
            reaction_children = [
                child
                for child in node.get("children", [])
                if child["type"] == "reaction"
            ]
            if len(reaction_children) > 1:
                print(f"Non-linear molecule (>1 reaction paths): {node['smiles']}")
                return False

        # Recursively check all children
        for child in node.get("children", []):
            if not is_linear_synthesis(child):
                return False

        return True

    is_linear = is_linear_synthesis(route)
    print(f"Synthesis is linear: {is_linear}")
    return is_linear
