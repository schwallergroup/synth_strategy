#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects early-stage introduction of cyclopropyl-containing groups.
    """
    cyclopropyl_introduced_early = False

    def dfs_traverse(node):
        nonlocal cyclopropyl_introduced_early

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an early-stage reaction (depth >= 5)
            if "depth" in node["metadata"] and node["metadata"].get("depth", -1) >= 5:
                # Try to detect cyclopropyl introduction
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    # SMARTS pattern for cyclopropyl group
                    cyclopropyl_pattern = Chem.MolFromSmarts("[CH]1[CH2][CH2]1")

                    # Check if product has cyclopropyl but reactants don't
                    if product_mol and product_mol.HasSubstructMatch(cyclopropyl_pattern):
                        print("Detected early-stage cyclopropyl introduction")
                        cyclopropyl_introduced_early = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cyclopropyl_introduced_early
