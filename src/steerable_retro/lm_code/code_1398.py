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
    Detects a strategy involving coupling of fragments where at least one contains
    a stereocenter that is preserved throughout the synthesis.
    """
    found_stereocenter = False
    found_coupling_with_stereocenter = False

    def has_stereocenter(smiles):
        """Check if molecule contains a stereocenter"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return "@" in smiles  # Simple check for @ symbol in SMILES
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_stereocenter, found_coupling_with_stereocenter

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for stereocenter in reactants and product
                stereo_reactant = False
                for r in reactants:
                    if has_stereocenter(r):
                        stereo_reactant = True
                        found_stereocenter = True

                product_has_stereo = has_stereocenter(product)

                # Check for coupling reaction (multiple reactants) where stereocenter is preserved
                if len(reactants) >= 2 and stereo_reactant and product_has_stereo:
                    print(f"Found stereoselective coupling at depth {depth}")
                    found_coupling_with_stereocenter = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if strategy criteria are met
    return found_stereocenter and found_coupling_with_stereocenter
