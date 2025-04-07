#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if a stereocenter is maintained throughout the entire synthesis.
    """
    # Track stereocenters by their atom mapping IDs
    stereocenters_by_depth = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol:
                # Get chiral atoms in this molecule
                chiral_atoms = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                if chiral_atoms:
                    if depth not in stereocenters_by_depth:
                        stereocenters_by_depth[depth] = set()
                    stereocenters_by_depth[depth].update(
                        [atom_idx for atom_idx, _ in chiral_atoms]
                    )
                    print(
                        f"Depth {depth}: Found {len(chiral_atoms)} stereocenters in {mol_smiles}"
                    )

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if rsmi:
                    # Extract product and reactants from atom-mapped reaction SMILES
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check stereocenters in product and map to reactants
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Get atom mapping from product
                        product_atom_maps = {}
                        for atom in product_mol.GetAtoms():
                            map_num = atom.GetAtomMapNum()
                            if map_num > 0:
                                product_atom_maps[atom.GetIdx()] = map_num

                        # Get chiral atoms in product
                        product_chiral_atoms = Chem.FindMolChiralCenters(
                            product_mol, includeUnassigned=False
                        )

                        # Track mapped stereocenters
                        for atom_idx, stereo in product_chiral_atoms:
                            if atom_idx in product_atom_maps:
                                map_num = product_atom_maps[atom_idx]
                                if depth not in stereocenters_by_depth:
                                    stereocenters_by_depth[depth] = set()
                                stereocenters_by_depth[depth].add(map_num)
                                print(
                                    f"Depth {depth}: Found stereocenter with map ID {map_num} in product"
                                )

                        # Check if stereocenters in reactants match those in product
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Get atom mapping from reactant
                                reactant_atom_maps = {}
                                for atom in reactant_mol.GetAtoms():
                                    map_num = atom.GetAtomMapNum()
                                    if map_num > 0:
                                        reactant_atom_maps[atom.GetIdx()] = map_num

                                # Get chiral atoms in reactant
                                reactant_chiral_atoms = Chem.FindMolChiralCenters(
                                    reactant_mol, includeUnassigned=False
                                )

                                # Track mapped stereocenters
                                for atom_idx, stereo in reactant_chiral_atoms:
                                    if atom_idx in reactant_atom_maps:
                                        map_num = reactant_atom_maps[atom_idx]
                                        if depth + 1 not in stereocenters_by_depth:
                                            stereocenters_by_depth[depth + 1] = set()
                                        stereocenters_by_depth[depth + 1].add(map_num)
                                        print(
                                            f"Depth {depth+1}: Found stereocenter with map ID {map_num} in reactant"
                                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if there are stereocenters at all depths
    if not stereocenters_by_depth:
        print("No stereocenters found in the synthesis route")
        return False

    # Find stereocenters that appear at all depths
    all_depths = set(range(max_depth + 1))
    depths_with_stereocenters = set(stereocenters_by_depth.keys())

    # Check if stereocenters exist at all depths
    if depths_with_stereocenters != all_depths:
        print(
            f"Stereocenters not present at all depths. Missing at depths: {all_depths - depths_with_stereocenters}"
        )
        return False

    # Check if at least one stereocenter exists at each depth
    has_maintained_stereocenter = all(
        len(stereocenters_by_depth.get(depth, [])) > 0 for depth in range(max_depth + 1)
    )

    print(
        f"Maintained stereocenter throughout synthesis: {has_maintained_stereocenter}"
    )
    if has_maintained_stereocenter:
        print(f"Stereocenters present at all depths")

    return has_maintained_stereocenter
