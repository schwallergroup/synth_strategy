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
    This function detects a linear synthesis with sequential heteroatom bond formations
    (C-O, C-N) connecting aromatic fragments.
    """
    # Track reaction count and linearity
    reaction_count = 0
    is_linear = True
    heteroatom_bond_formations = 0

    def dfs_traverse(node, depth=0, parent_children_count=1):
        nonlocal reaction_count, is_linear, heteroatom_bond_formations

        if node["type"] == "reaction":
            reaction_count += 1

            # Check if this is a branching point (convergent synthesis)
            if len(node.get("children", [])) > 1:
                # More than one reactant means potentially convergent
                if parent_children_count > 1:
                    is_linear = False

            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants = parts[0].split(".")
                    product = parts[2]

                    # Check for C-O bond formation
                    if "[O:" in rsmi and (
                        ("[CH2:" in rsmi and "[O:" in rsmi) or ("[OH:" in rsmi and "[c:" in rsmi)
                    ):
                        print(f"C-O bond formation detected at depth {depth}")
                        heteroatom_bond_formations += 1

                    # Check for C-N bond formation
                    if "[NH:" in rsmi and "[c:" in rsmi:
                        print(f"C-N bond formation detected at depth {depth}")
                        heteroatom_bond_formations += 1

        # Process children
        children_count = len(node.get("children", []))
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, children_count)

    # Start traversal
    dfs_traverse(route)

    # Criteria for this strategy
    result = reaction_count >= 2 and is_linear and heteroatom_bond_formations >= 2

    print(f"Linear synthesis with heteroatom bonds: {result}")
    print(f"  - Reaction count: {reaction_count}")
    print(f"  - Is linear: {is_linear}")
    print(f"  - Heteroatom bond formations: {heteroatom_bond_formations}")

    return result
