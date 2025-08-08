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
    This function detects convergent synthesis approach with multiple fragment couplings.
    """
    fragment_coupling_count = 0

    def dfs_traverse(node):
        nonlocal fragment_coupling_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If reaction has multiple reactants, it's potentially a fragment coupling
                if len(reactants) >= 2:
                    # Check if reactants are substantial fragments (not just reagents)
                    substantial_fragments = 0
                    for r in reactants:
                        mol = Chem.MolFromSmiles(r)
                        if (
                            mol and mol.GetNumHeavyAtoms() > 6
                        ):  # Arbitrary threshold for "substantial"
                            substantial_fragments += 1

                    if substantial_fragments >= 2:
                        fragment_coupling_count += 1
                        print(f"Detected fragment coupling (count: {fragment_coupling_count})")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Return True if multiple fragment couplings are detected
    return fragment_coupling_count >= 2
