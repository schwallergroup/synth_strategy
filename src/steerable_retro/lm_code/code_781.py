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
    Detects late-stage formation of sulfonamide group (in first half of synthesis).
    """
    sulfonamide_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains sulfonamide but reactants don't
            sulfonamide_pattern = Chem.MolFromSmarts("[NH]-[S](=[O])(=[O])-[c]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                # Check if reactants don't have sulfonamide
                has_sulfonamide_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(sulfonamide_pattern):
                        has_sulfonamide_in_reactants = True
                        break

                if not has_sulfonamide_in_reactants:
                    sulfonamide_formation_depth = depth
                    print(f"Sulfonamide formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if sulfonamide formation occurs in first half of synthesis (late stage)
    if sulfonamide_formation_depth is not None and sulfonamide_formation_depth <= max_depth / 2:
        return True
    return False
