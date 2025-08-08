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
    This function detects if the synthesis follows a convergent approach
    where two or more complex fragments are combined in the final step.
    """
    is_convergent = False

    def dfs_traverse(node):
        nonlocal is_convergent

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Get depth information
            depth = 0
            if "depth" in node:
                depth = node["depth"]
            elif "metadata" in node and "depth" in node["metadata"]:
                depth = node["metadata"]["depth"]

            # Check if this is the final reaction (depth 0)
            if depth == 0:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # If there are 2 or more complex reactants in the final step
                if len(reactants) >= 2:
                    # Check if reactants are complex (not simple reagents)
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > 6:  # Arbitrary threshold for "complex"
                            complex_reactants += 1

                    if complex_reactants >= 2:
                        print(
                            "Detected convergent synthesis with multiple complex fragments in final step"
                        )
                        is_convergent = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return is_convergent
