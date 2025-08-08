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
    This function detects if both bromide and fluoride substituents are maintained
    throughout the synthesis (from starting materials to final product).
    """
    # Find starting materials and final product
    starting_materials = []
    final_product = None

    # Track the depth of each node to identify the final product (lowest depth)
    def find_nodes(node, depth=0):
        nonlocal final_product

        if node["type"] == "mol":
            # If it's a starting material (in_stock), add to starting materials
            if node.get("in_stock", False):
                starting_materials.append(node)

            # Update final product if this is at a lower depth (or first molecule found)
            if final_product is None or depth < final_product["depth"]:
                final_product = {"smiles": node["smiles"], "depth": depth}

        # Traverse children with increased depth
        for child in node.get("children", []):
            find_nodes(child, depth + 1)

    find_nodes(route)

    # If we couldn't identify starting materials or final product, return False
    if not starting_materials or final_product is None:
        print("Could not identify starting materials or final product")
        return False

    # Check if both bromide and fluoride are present in any starting material
    bromide_in_starting = False
    fluoride_in_starting = False

    for material in starting_materials:
        # Check for bromide
        if "Br" in material["smiles"]:
            bromide_in_starting = True

        # Check for fluoride
        if "F" in material["smiles"]:
            fluoride_in_starting = True

    print(f"Bromide in starting materials: {bromide_in_starting}")
    print(f"Fluoride in starting materials: {fluoride_in_starting}")

    # If neither halogen is in starting materials, return True (nothing to maintain)
    if not bromide_in_starting and not fluoride_in_starting:
        print("No halogens found in starting materials")
        return True

    # Check if the halogens from starting materials are present in final product
    final_smiles = final_product["smiles"]
    bromide_in_final = "Br" in final_smiles
    fluoride_in_final = "F" in final_smiles

    print(f"Bromide in final product: {bromide_in_final}")
    print(f"Fluoride in final product: {fluoride_in_final}")

    # Halogens are maintained if they're in both starting materials and final product
    bromide_maintained = (not bromide_in_starting) or bromide_in_final
    fluoride_maintained = (not fluoride_in_starting) or fluoride_in_final

    if not bromide_maintained:
        print("Bromide was present in starting materials but lost in final product")
    if not fluoride_maintained:
        print("Fluoride was present in starting materials but lost in final product")

    return bromide_maintained and fluoride_maintained
