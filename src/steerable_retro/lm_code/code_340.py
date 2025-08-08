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
    This function detects convergent synthesis where three or more fragments
    are combined in the final step.
    """
    found_convergent = False

    # In retrosynthetic analysis, the final step is the first reaction node
    # we encounter when starting from the root
    def find_final_step(node):
        if node["type"] == "reaction":
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")

                    # Check if it combines 3 or more fragments
                    if len(reactants) >= 3:
                        print(
                            f"Found convergent synthesis with {len(reactants)} fragments in final step"
                        )
                        return True
            except Exception as e:
                print(f"Error processing reaction node: {e}")
        return False

    # If the root is a molecule node, check its immediate reaction children
    if route["type"] == "mol":
        for child in route.get("children", []):
            if find_final_step(child):
                found_convergent = True
                break
    # If the root is already a reaction node, check it directly
    elif route["type"] == "reaction":
        found_convergent = find_final_step(route)

    return found_convergent
