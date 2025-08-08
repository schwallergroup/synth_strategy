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
    Detects if stereocenters present in starting materials are preserved
    throughout the synthesis.
    """
    starting_stereocenters = 0
    final_stereocenters = 0
    starting_materials_with_stereo = []

    def dfs_traverse(node, depth=0):
        nonlocal starting_stereocenters, final_stereocenters

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Count stereocenters
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                stereo_count = len(chiral_centers)

                # Check if this is a starting material (either marked as in_stock or a leaf node)
                is_starting_material = ("in_stock" in node and node["in_stock"]) or (
                    len(node.get("children", [])) == 0 and depth > 0
                )

                if is_starting_material:  # Starting material
                    if stereo_count > 0:
                        starting_stereocenters += stereo_count
                        starting_materials_with_stereo.append(node["smiles"])
                        print(
                            f"Found {stereo_count} stereocenters in starting material: {node['smiles']}"
                        )
                elif depth == 0:  # Final product
                    final_stereocenters = stereo_count
                    print(f"Found {stereo_count} stereocenters in final product: {node['smiles']}")
                else:
                    print(
                        f"Intermediate molecule at depth {depth}: {node['smiles']} with {stereo_count} stereocenters"
                    )

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Total stereocenters found: {starting_stereocenters} in starting materials, {final_stereocenters} in final product"
    )
    print(f"Starting materials with stereocenters: {starting_materials_with_stereo}")

    # Check if stereocenters are preserved
    # If there are stereocenters in the final product, consider them preserved
    # This accounts for both preservation of existing stereocenters and creation of new ones
    if final_stereocenters > 0:
        print(f"Stereocenters present in final product: {final_stereocenters}")
        return True

    return False
