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
    Detects if the route preserves a stereocenter from starting materials to final product.
    """
    # Track stereocenters by atom mapping
    stereocenter_maps = {}
    # Track depths where each mapped stereocenter appears
    stereocenter_depths = {}

    def dfs_traverse(node, depth=0):
        # Set the depth for the current node
        node["depth"] = depth

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            try:
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Find stereocenters in the molecule
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)

                    if chiral_centers:
                        print(
                            f"Found {len(chiral_centers)} stereocenters at depth {depth} in {mol_smiles}"
                        )

                        # Get atom mapping from parent reaction if not the root node
                        atom_maps = {}
                        if depth > 0 and "children" in route:
                            # Find the reaction node that produced this molecule
                            for reaction_node in route.get("children", []):
                                if (
                                    reaction_node["type"] == "reaction"
                                    and "metadata" in reaction_node
                                ):
                                    try:
                                        rsmi = reaction_node["metadata"].get("rsmi", "")
                                        if rsmi:
                                            product = rsmi.split(">")[-1]
                                            if Chem.MolFromSmiles(product) and Chem.MolToSmiles(
                                                Chem.MolFromSmiles(product)
                                            ) == Chem.MolToSmiles(mol):
                                                # Extract atom mapping from the product
                                                product_mol = Chem.MolFromSmiles(product)
                                                for atom in product_mol.GetAtoms():
                                                    if atom.GetAtomMapNum() > 0:
                                                        atom_maps[atom.GetIdx()] = (
                                                            atom.GetAtomMapNum()
                                                        )
                                    except Exception as e:
                                        print(f"Error processing reaction SMILES: {e}")

                        # Record stereocenters with their atom mapping
                        for atom_idx, stereo in chiral_centers:
                            # If we have atom mapping, use it to track the stereocenter
                            if atom_idx in atom_maps:
                                map_num = atom_maps[atom_idx]
                                if map_num not in stereocenter_maps:
                                    stereocenter_maps[map_num] = []
                                stereocenter_maps[map_num].append((depth, stereo))

                                if map_num not in stereocenter_depths:
                                    stereocenter_depths[map_num] = set()
                                stereocenter_depths[map_num].add(depth)
                            else:
                                # For the root node or if no mapping is available
                                if depth not in stereocenter_maps:
                                    stereocenter_maps[depth] = []
                                stereocenter_maps[depth].append((atom_idx, stereo))
            except Exception as e:
                print(f"Error processing molecule SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Stereocenter maps: {stereocenter_maps}")
    print(f"Stereocenter depths: {stereocenter_depths}")

    # Check if any stereocenter appears at both early and late stages
    for map_num, depths in stereocenter_depths.items():
        depths_list = list(depths)
        if len(depths_list) >= 2 and min(depths_list) <= 1 and max(depths_list) >= 2:
            print(
                f"Stereocenter with map {map_num} preserved from depth {min(depths_list)} to {max(depths_list)}"
            )
            return True

    # Fallback to simpler check if atom mapping approach doesn't find preserved stereocenters
    stereocenters_by_depth = {}

    def count_stereocenters(node):
        if node["type"] == "mol" and "smiles" in node:
            depth = node.get("depth", None)
            if depth is not None:
                smiles = node["smiles"]
                # Check for stereochemistry markers in SMILES
                if "@" in smiles:
                    if depth not in stereocenters_by_depth:
                        stereocenters_by_depth[depth] = 0
                    stereocenters_by_depth[depth] += smiles.count("@")
                    print(f"Counted {smiles.count('@')} stereocenters at depth {depth}")

        for child in node.get("children", []):
            count_stereocenters(child)

    count_stereocenters(route)

    print(f"Stereocenters by depth (fallback): {stereocenters_by_depth}")

    # If we have stereocenters at both early and late stages
    depths = list(stereocenters_by_depth.keys())
    if len(depths) >= 2 and min(depths) <= 1 and max(depths) >= 2:
        print(
            f"Found stereocenters at both early (depth {max(depths)}) and late (depth {min(depths)}) stages"
        )
        return True

    return False
