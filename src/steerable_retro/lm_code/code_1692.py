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
    This function detects if an ester is converted to an amide during the synthesis.
    """
    ester_to_amide_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_to_amide_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester in reactants
                ester_pattern = Chem.MolFromSmarts("C(=O)OC")
                ester_found = False
                for reactant in reactants:
                    if (
                        reactant
                        and Chem.MolFromSmiles(reactant)
                        and Chem.MolFromSmiles(reactant).HasSubstructMatch(ester_pattern)
                    ):
                        ester_found = True
                        break

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("C(=O)N")
                amide_found = False
                if (
                    product
                    and Chem.MolFromSmiles(product)
                    and Chem.MolFromSmiles(product).HasSubstructMatch(amide_pattern)
                ):
                    amide_found = True

                if ester_found and amide_found:
                    print("Ester to amide conversion detected at depth", depth)
                    ester_to_amide_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return ester_to_amide_detected
