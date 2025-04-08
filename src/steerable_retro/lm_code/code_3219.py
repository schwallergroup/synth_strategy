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
    This function detects aromatic iodination in the synthetic route.
    """
    iodination_detected = False

    def dfs_traverse(node):
        nonlocal iodination_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Try to detect iodination pattern
            try:
                p_mol = Chem.MolFromSmiles(product)
                aromatic_iodine_pattern = Chem.MolFromSmarts("[c]-[I]")

                # Check if product has aromatic iodine
                if p_mol and p_mol.HasSubstructMatch(aromatic_iodine_pattern):
                    # Check if any reactant contains iodine source
                    for r in reactants:
                        if "I" in r or "ClI" in r:  # Simple check for iodine-containing reagents
                            print("Aromatic iodination detected in reaction:", rsmi)
                            iodination_detected = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return iodination_detected
