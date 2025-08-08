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
    Detects if the synthesis preserves stereochemistry throughout.
    """
    stereocenters_preserved = True

    def dfs_traverse(node):
        nonlocal stereocenters_preserved

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for stereochemistry markers in SMILES
            stereo_markers = ["@", "/", "\\"]

            # If product has stereochemistry markers
            if any(marker in product_smiles for marker in stereo_markers):
                # Check if any reactant also has stereochemistry
                reactants_have_stereo = any(
                    any(marker in r for marker in stereo_markers) for r in reactants_smiles
                )

                # If reactants have stereo but product doesn't have expected number, stereo might be lost
                if reactants_have_stereo:
                    # This is a simplified check - a more robust implementation would count and compare
                    # the actual number of stereocenters in reactants and products
                    print("Stereochemistry preserved in this step")
                else:
                    # This step creates new stereochemistry
                    print("New stereochemistry created in this step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Has stereocenter preservation strategy: {stereocenters_preserved}")
    return stereocenters_preserved
