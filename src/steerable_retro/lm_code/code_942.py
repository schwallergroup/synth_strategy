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
    This function detects Wittig-type reactions for C=C bond formation early in the synthesis.
    """
    wittig_detected = False
    max_depth = 0
    wittig_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal wittig_detected, max_depth, wittig_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phosphorus-containing reactant (typical for Wittig)
                phosphorus_present = False
                for reactant in reactants:
                    if "P" in reactant:
                        phosphorus_present = True
                        break

                # Check if product contains C=C bond
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
                    if phosphorus_present:
                        wittig_detected = True
                        if wittig_depth is None or depth > wittig_depth:
                            wittig_depth = depth
                            print(f"Wittig-type reaction detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if Wittig reaction occurs in the second half of synthesis (early in forward direction)
    if wittig_depth is not None and wittig_depth > max_depth / 2:
        return True
    return False
