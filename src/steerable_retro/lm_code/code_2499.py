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
    Detects if stereochemistry is preserved throughout the synthesis.

    This function tracks stereocenters through the synthesis route using atom mapping
    in reaction SMILES. It verifies that existing stereocenters maintain their
    configuration throughout the synthesis process.
    """
    # Track stereocenters through the synthesis
    stereo_tracking = {}  # Maps (depth, atom_map) to stereochemistry

    def get_chiral_atoms_with_mapping(mol):
        """Extract chiral atoms with their mapping IDs"""
        chiral_atoms = {}
        if not mol:
            return chiral_atoms

        # Find all chiral centers
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        # Get atom mapping for each chiral center
        for atom_idx, chirality in chiral_centers:
            atom = mol.GetAtomWithIdx(atom_idx)
            map_num = atom.GetAtomMapNum()

            if map_num > 0:  # Valid atom mapping
                # Store the chiral tag (R or S configuration)
                chiral_tag = atom.GetChiralTag()
                chiral_atoms[map_num] = (atom_idx, chiral_tag)
                print(
                    f"Found chiral atom with map {map_num}, chirality: {chirality}, tag: {chiral_tag}"
                )

        return chiral_atoms

    def dfs_traverse(node, depth=0):
        nonlocal stereo_tracking

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"], sanitize=True)
                if mol:
                    # Get chiral atoms with their mapping
                    chiral_atoms = get_chiral_atoms_with_mapping(mol)

                    # Store stereochemistry information at this depth
                    for map_num, (atom_idx, chiral_tag) in chiral_atoms.items():
                        stereo_tracking[(depth, map_num)] = chiral_tag

                    print(
                        f"Depth {depth}, SMILES: {node['smiles']}, Chiral atoms: {len(chiral_atoms)}"
                    )
            except Exception as e:
                print(f"Error processing SMILES in stereochemistry detection: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                # Extract reactants and products with atom mapping
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                # Process reactants to find stereocenters
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                if reactants_mol:
                    reactant_chiral_atoms = get_chiral_atoms_with_mapping(reactants_mol)
                    print(
                        f"Depth {depth}, Reaction, Reactant chiral atoms: {len(reactant_chiral_atoms)}"
                    )

                # Process products to find stereocenters
                products_mol = Chem.MolFromSmiles(products_smiles)
                if products_mol:
                    product_chiral_atoms = get_chiral_atoms_with_mapping(products_mol)
                    print(
                        f"Depth {depth}, Reaction, Product chiral atoms: {len(product_chiral_atoms)}"
                    )

                    # Check if stereocenters are preserved in this reaction
                    for map_num in reactant_chiral_atoms:
                        if map_num in product_chiral_atoms:
                            # If the atom is still chiral in the product, check if configuration is preserved
                            reactant_chiral_tag = reactant_chiral_atoms[map_num][1]
                            product_chiral_tag = product_chiral_atoms[map_num][1]

                            # Store product stereochemistry
                            stereo_tracking[(depth + 1, map_num)] = product_chiral_tag

                            # Check if stereochemistry changed
                            if (
                                reactant_chiral_tag != product_chiral_tag
                                and reactant_chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED
                                and product_chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED
                            ):
                                print(
                                    f"Stereochemistry changed for atom map {map_num} at depth {depth}"
                                )
                        else:
                            # Stereocenter was lost
                            print(f"Stereocenter with map {map_num} was lost at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if any stereocenters exist in the synthesis route
    all_stereocenters = stereo_tracking.keys()

    # Check if stereocenters are preserved throughout the synthesis
    preserved = True

    # Get all unique atom mappings that appear as stereocenters
    all_atom_maps = set(map_num for _, map_num in all_stereocenters)

    for map_num in all_atom_maps:
        # Get all depths where this atom map appears
        depths = sorted([d for d, m in stereo_tracking.keys() if m == map_num])

        if len(depths) > 1:  # This stereocenter appears in multiple steps
            # Check if stereochemistry is consistent
            first_tag = None
            for d in depths:
                if (d, map_num) in stereo_tracking:
                    current_tag = stereo_tracking[(d, map_num)]
                    if first_tag is None:
                        first_tag = current_tag
                    elif (
                        current_tag != first_tag
                        and current_tag != Chem.ChiralType.CHI_UNSPECIFIED
                        and first_tag != Chem.ChiralType.CHI_UNSPECIFIED
                    ):
                        preserved = False
                        print(
                            f"Stereochemistry not preserved for atom map {map_num} across depths {depths}"
                        )
                        break

    # Final check: if we have stereocenters and they're preserved
    if preserved and all_atom_maps:
        print("Stereochemistry is preserved throughout synthesis")
        return True
    elif not all_atom_maps:
        print("No stereocenters found in the synthesis route")
        return False
    else:
        print("Stereochemistry is not preserved")
        return False
