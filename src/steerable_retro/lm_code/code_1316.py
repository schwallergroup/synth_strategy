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
    This function detects the use of tosylate as a leaving group in the synthesis.
    """
    tosylate_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6]([S](=[O])(=[O])[O])[#6][#6]1")

    strategy_detected = False

    def dfs_traverse(node):
        nonlocal strategy_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check if any reactant contains tosylate group
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(tosylate_pattern):
                        print("Detected tosylate leaving group strategy")
                        strategy_detected = True
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return strategy_detected
