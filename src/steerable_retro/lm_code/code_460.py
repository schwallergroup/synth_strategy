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
    Detects if the synthesis maintains stereochemistry throughout the route.

    Returns True if all stereocenters in the target molecule are preserved
    throughout the synthesis route.
    """
    # Track stereocenters by depth and atom mapping
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        # Handle molecule nodes (including target molecule)
        if node["type"] == "mol":
            mol_smiles = node.get("smiles", "")
            if not mol_smiles:
                return

            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                print(f"Warning: Could not parse molecule SMILES: {mol_smiles}")
                return

            # Find chiral centers in the molecule
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
            if chiral_centers:
                stereocenters_by_depth[depth] = {
                    "count": len(chiral_centers),
                    "centers": chiral_centers,
                    "smiles": mol_smiles,
                }
                print(
                    f"Depth {depth} (molecule): Found {len(chiral_centers)} stereocenters in {mol_smiles}"
                )

        # Handle reaction nodes
        elif node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            product_smiles = rsmi.split(">")[-1]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not product:
                print(f"Warning: Could not parse product SMILES: {product_smiles}")
                return

            # Find chiral centers in the product
            chiral_centers = Chem.FindMolChiralCenters(product, includeUnassigned=False)
            if chiral_centers:
                stereocenters_by_depth[depth] = {
                    "count": len(chiral_centers),
                    "centers": chiral_centers,
                    "smiles": product_smiles,
                }
                print(
                    f"Depth {depth} (reaction product): Found {len(chiral_centers)} stereocenters in {product_smiles}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If no stereocenters found, return False
    if not stereocenters_by_depth:
        print("No stereocenters found in the synthesis route")
        return False

    # Get stereocenters in the target molecule (depth 0)
    target_stereocenters = stereocenters_by_depth.get(0, None)

    # If target has no stereocenters but intermediates do, that's still a valid strategy
    # (stereocenters are introduced and then maintained)
    if not target_stereocenters:
        min_depth = min(stereocenters_by_depth.keys())
        target_stereocenters = stereocenters_by_depth[min_depth]
        print(f"Target molecule has no stereocenters, using depth {min_depth} as reference")

    target_count = target_stereocenters["count"]

    # Check if the number of stereocenters is maintained or increases throughout the synthesis
    # (we allow for introducing new stereocenters in later stages)
    maintained = True
    for depth, info in sorted(stereocenters_by_depth.items()):
        if depth > 0 and info["count"] < target_count:
            print(
                f"Stereochemistry not maintained at depth {depth}: {info['count']} < {target_count}"
            )
            maintained = False
            break

    if maintained:
        print("Stereochemistry maintained throughout synthesis")

    return maintained
