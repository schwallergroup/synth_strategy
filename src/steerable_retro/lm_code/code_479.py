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
    Detects if the synthesis follows a late-stage convergent approach where
    two major fragments are combined in the final step.
    """
    final_coupling_detected = False

    def dfs_traverse(node):
        nonlocal final_coupling_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this is the final reaction (depth 0)
            if node.get("depth", 0) == 0:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If we have multiple reactants at the final step, it's convergent
                if len(reactants) >= 2:
                    # Check if reactants are substantial fragments (not just small reagents)
                    substantial_fragments = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > 10:  # Arbitrary threshold for "substantial"
                            substantial_fragments += 1

                    if substantial_fragments >= 2:
                        print(
                            "Late-stage convergent synthesis detected: multiple substantial fragments combined in final step"
                        )
                        final_coupling_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return final_coupling_detected
