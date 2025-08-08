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
    Detects if the synthesis preserves stereochemistry throughout the route.
    """
    has_stereocenter = False
    preserves_stereocenter = True

    def dfs_traverse(node, depth=0):
        nonlocal has_stereocenter, preserves_stereocenter

        if node["type"] == "mol":
            smiles = node["smiles"]
            # Check if molecule has stereocenter
            if "@" in smiles:
                has_stereocenter = True
                print(f"Stereocenter detected at depth {depth}: {smiles}")

        elif node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reaction preserves stereochemistry
            if "@" in product:
                has_stereocenter = True
                # Check if all reactants with @ maintain their stereochemistry
                for reactant in reactants:
                    if "@" in reactant and not any(stereo in product for stereo in ["@H", "@@H"]):
                        preserves_stereocenter = False
                        print(f"Stereocenter not preserved in reaction at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_stereocenter and preserves_stereocenter
