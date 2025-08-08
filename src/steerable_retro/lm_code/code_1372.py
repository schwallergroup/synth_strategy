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


def main(route, debug=False):
    """
    Detects if the synthesis follows a linear strategy (no convergent steps with multiple complex fragments).

    A linear synthesis typically builds a molecule by sequential addition of small fragments to a growing core.
    Convergent synthesis involves separate complex fragments being joined in late-stage reactions.

    Args:
        route: The synthesis route tree
        debug: Whether to print debug information

    Returns:
        bool: True if the synthesis is linear, False if convergent steps are detected
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Handle empty reactants
            if not reactants_part:
                return

            reactants = reactants_part.split(".")

            # Skip if only one reactant (definitely linear)
            if len(reactants) <= 1:
                return

            # Analyze reactant complexity
            reactant_complexities = []
            total_atoms = 0

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        num_atoms = mol.GetNumAtoms()
                        num_rings = len(mol.GetSSSR())

                        # Calculate complexity score based on atoms and rings
                        complexity = num_atoms + (num_rings * 2)
                        reactant_complexities.append(complexity)
                        total_atoms += num_atoms
                except Exception as e:
                    if debug:
                        print(f"Error processing reactant {reactant}: {e}")
                    continue

            # Sort complexities in descending order
            reactant_complexities.sort(reverse=True)

            # Check for convergent pattern - multiple significant fragments
            if len(reactant_complexities) >= 2:
                # If the second most complex reactant is significant compared to the total
                # (more than 25% of total atoms or complexity > 12)
                if len(reactant_complexities) >= 2 and (
                    (reactant_complexities[1] > 12)
                    or (total_atoms > 0 and reactant_complexities[1] / total_atoms > 0.25)
                ):
                    is_linear = False
                    if debug:
                        print(
                            f"Found convergent step with reactant complexities: {reactant_complexities}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
