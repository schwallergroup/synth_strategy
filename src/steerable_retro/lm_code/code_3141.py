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
    Detects if the route maintains a stereocenter throughout the synthesis.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        nonlocal stereocenters_by_depth

        if node["type"] == "reaction":
            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for stereocenters in the product
            if "@" in product:  # Quick check for stereochemistry
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Count stereocenters
                    chiral_centers = Chem.FindMolChiralCenters(product_mol, includeUnassigned=False)
                    if chiral_centers:
                        stereocenters_by_depth[depth] = len(chiral_centers)
                        print(f"Found {len(chiral_centers)} stereocenters at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if stereocenters are maintained across multiple steps
    depths_with_stereocenters = list(stereocenters_by_depth.keys())
    result = len(depths_with_stereocenters) >= 2

    print(f"Stereocenter preservation strategy detected: {result}")
    return result
