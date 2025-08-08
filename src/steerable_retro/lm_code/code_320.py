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
    This function detects a strategy involving coupling of aromatic fragments.
    """
    aromatic_couplings = 0

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_couplings

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants = parts[0].split(".")
                    product = parts[2]

                    # Check if reaction involves aromatic systems
                    aromatic_count_reactants = sum(1 for r in reactants if "[c:" in r or "c1" in r)
                    aromatic_in_product = "[c:" in product or "c1" in product

                    if aromatic_count_reactants >= 1 and aromatic_in_product:
                        # Check for bond formation between aromatic systems
                        if ("[c:" in rsmi and "[O:" in rsmi) or ("[c:" in rsmi and "[N:" in rsmi):
                            aromatic_couplings += 1
                            print(f"Aromatic fragment coupling detected at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Criteria for this strategy
    result = aromatic_couplings >= 2

    print(f"Aromatic fragment coupling strategy: {result}")
    print(f"  - Aromatic couplings: {aromatic_couplings}")

    return result
