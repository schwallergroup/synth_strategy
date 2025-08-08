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
    Detects if the synthesis involves sequential functionalization of an existing scaffold.
    """
    # Count different types of functionalizations
    functionalization_types = {
        "protection": False,
        "oxidation": False,
        "c_c_bond_formation": False,
        "functional_group_interconversion": False,
    }

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for protection (OH to OCH3)
            if "[OH:" in "".join(reactants) and "[O:" in product and "[CH3:" in product:
                functionalization_types["protection"] = True

            # Check for oxidation (S to SO2)
            if "[S:" in "".join(reactants) and "[S:" in product and "=[O:" in product:
                functionalization_types["oxidation"] = True

            # Check for C-C bond formation
            if "[c:" in product and "[C:" in product and "=[O:" in product:
                functionalization_types["c_c_bond_formation"] = True

            # Check for functional group interconversion (ketone to acid)
            if (
                "[C:" in product
                and "=[O:" in product
                and "[OH:" in product
                and "[CH3:" in "".join(reactants)
            ):
                functionalization_types["functional_group_interconversion"] = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Count how many different types of functionalizations are present
    functionalization_count = sum(1 for value in functionalization_types.values() if value)

    # If at least 3 different types of functionalizations, consider it sequential functionalization
    if functionalization_count >= 3:
        print(
            f"Sequential functionalization detected with {functionalization_count} different types"
        )
        return True
    return False
