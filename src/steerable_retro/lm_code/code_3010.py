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
    This function detects a linear synthesis strategy where a CF3-containing building block
    is introduced in the late stage of the synthesis.
    """
    cf3_introduction_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal cf3_introduction_depth, max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for CF3 introduction
            cf3_pattern = Chem.MolFromSmarts("[#6]-[C]([F])([F])[F]")
            product_mol = Chem.MolFromSmiles(product)

            # Check if product has CF3
            if product_mol and product_mol.HasSubstructMatch(cf3_pattern):
                # Check if any reactant doesn't have CF3
                cf3_in_all_reactants = True
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and not reactant_mol.HasSubstructMatch(cf3_pattern):
                        cf3_in_all_reactants = False
                        break

                if not cf3_in_all_reactants:
                    print(f"Found CF3 introduction at depth {depth}")
                    cf3_introduction_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if CF3 was introduced in the late stage (first half of synthesis)
    if cf3_introduction_depth is not None and cf3_introduction_depth <= max_depth / 2:
        print(f"CF3 introduced in late stage (depth {cf3_introduction_depth} of {max_depth})")
        return True
    return False
