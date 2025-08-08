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
    This function detects a linear synthesis strategy (as opposed to convergent).

    A linear synthesis strategy builds the target molecule by adding one significant
    building block at a time in a sequential manner. In contrast, a convergent strategy
    combines multiple complex fragments in later stages.
    """
    # Track if we have any convergent steps (reactions with multiple complex reactants)
    has_convergent_steps = False

    def dfs_traverse(node):
        nonlocal has_convergent_steps

        if node["type"] == "reaction" and "children" in node:
            # Count significant non-in-stock reactants for this reaction
            significant_reactants = 0

            for child in node["children"]:
                if child["type"] == "mol" and not child.get("in_stock", True):
                    # Parse the molecule
                    mol = Chem.MolFromSmiles(child["smiles"])
                    if mol is None:
                        print(f"Warning: Could not parse SMILES: {child['smiles']}")
                        continue

                    # Consider a molecule significant if it has more than 7 atoms
                    # and is not marked as in_stock
                    atom_count = len(mol.GetAtoms())
                    if atom_count > 7:
                        significant_reactants += 1

            # If a reaction has more than one significant reactant, it's convergent
            if significant_reactants > 1:
                has_convergent_steps = True
                print(f"Found convergent step with {significant_reactants} significant reactants")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # A synthesis is linear if it has no convergent steps
    is_linear = not has_convergent_steps
    print(f"Linear synthesis: {is_linear}")

    return is_linear
