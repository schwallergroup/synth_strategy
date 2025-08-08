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
    This function detects a linear synthesis strategy where each reaction
    has only one non-commercial reactant that leads to further reactions
    (indicating sequential addition).

    In a retrosynthetic tree, a linear synthesis means that at each step,
    only one of the reactants is derived from a previous reaction.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # Skip further processing if we already know it's not linear
        if not is_linear:
            return

        if node["type"] == "reaction":
            print(f"Examining reaction at depth {depth}")

            # Count non-commercial reactants that lead to further reactions
            complex_reactants = []

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    # Check if this non-commercial reactant leads to further reactions
                    has_reaction_children = False
                    for grandchild in child.get("children", []):
                        if grandchild["type"] == "reaction":
                            has_reaction_children = True
                            break

                    if has_reaction_children:
                        complex_reactants.append(child["smiles"])

            # If more than one complex reactant leads to further reactions, it's not linear
            if len(complex_reactants) > 1:
                is_linear = False
                print(
                    f"Non-linear reaction found at depth {depth} with {len(complex_reactants)} complex reactants"
                )
                for i, smiles in enumerate(complex_reactants):
                    print(f"  Complex reactant {i+1}: {smiles[:30]}...")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if is_linear:
        print("The synthesis route follows a linear strategy")
    else:
        print("The synthesis route does not follow a linear strategy")

    return is_linear
