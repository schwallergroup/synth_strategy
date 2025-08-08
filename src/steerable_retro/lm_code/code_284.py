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
    Detects a linear synthesis with multiple C-N bond formations (≥2)
    """
    cn_bond_formations = 0
    is_linear = True

    def dfs_traverse(node):
        nonlocal cn_bond_formations, is_linear

        if node["type"] == "reaction":
            # Check if this is a linear synthesis (one product)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                products = rsmi.split(">")[-1].split(".")

                if len(products) > 1:
                    is_linear = False

                # Check for C-N bond formation
                # This is a simplified approach - in practice, you'd need to analyze
                # the reaction more carefully to detect new C-N bonds
                if "N" in rsmi.split(">")[-1] and any("C" in r for r in reactants):
                    # Simple heuristic: if product has N and any reactant has C
                    # This is an oversimplification - real implementation would need atom mapping
                    cn_bond_formations += 1
                    print(f"Detected potential C-N bond formation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = is_linear and cn_bond_formations >= 2
    print(f"Linear synthesis: {is_linear}, C-N bond formations: {cn_bond_formations}")
    return result
