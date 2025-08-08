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
    This function detects if specific functional groups (halogens, ethers)
    are preserved throughout the synthesis.
    """
    # Track functional groups in final product
    final_product_groups = {"fluorine": False, "chlorine": False, "phenoxy": False}

    # Track if these groups are also in starting materials
    starting_material_groups = {"fluorine": False, "chlorine": False, "phenoxy": False}

    def check_functional_groups(smiles):
        """Check for presence of specific functional groups"""
        groups = {"fluorine": False, "chlorine": False, "phenoxy": False}

        # Check for fluorine
        if re.search(r"F|\[F\]", smiles):
            groups["fluorine"] = True

        # Check for chlorine
        if re.search(r"Cl|\[Cl\]", smiles):
            groups["chlorine"] = True

        # Check for phenoxy group (simplified check)
        if re.search(r"c[O]c|c\[O\]c", smiles):
            groups["phenoxy"] = True

        return groups

    def dfs_traverse(node, is_final=True):
        nonlocal final_product_groups, starting_material_groups

        if node["type"] == "mol":
            smiles = node["smiles"]
            groups = check_functional_groups(smiles)

            # If this is the final product, record its functional groups
            if is_final:
                final_product_groups = groups

            # If this is a starting material (in_stock), record its functional groups
            if node.get("in_stock", False):
                for group in groups:
                    if groups[group]:
                        starting_material_groups[group] = True

        # Traverse children (not final product anymore)
        for child in node.get("children", []):
            dfs_traverse(child, False)

    # Start traversal
    dfs_traverse(route)

    # Check which functional groups are preserved (present in both final product and starting materials)
    preserved = []
    for group in final_product_groups:
        if final_product_groups[group] and starting_material_groups[group]:
            preserved.append(group)
            print(f"Preserved functional group: {group}")

    # Return True if at least one functional group is preserved
    return len(preserved) > 0
