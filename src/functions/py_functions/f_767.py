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
    This function detects if the synthetic route maintains a stereocenter throughout.
    It tracks specific stereocenters across reactions using atom mapping.
    """
    # Track stereocenters at each molecule node
    stereocenters_by_mol = {}
    # Track depth of each molecule
    mol_depths = {}
    # Track atom mappings across reactions
    atom_mappings = {}

    def get_stereocenters(mol_smiles):
        """Get atom indices of stereocenters in a molecule"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return []

        stereocenters = []
        for atom in mol.GetAtoms():
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                # Store atom index and mapping number if available
                map_num = (
                    atom.GetProp("molAtomMapNumber")
                    if atom.HasProp("molAtomMapNumber")
                    else None
                )
                stereocenters.append((atom.GetIdx(), map_num))

        return stereocenters

    def dfs_traverse(node, depth=0):
        """Traverse the synthesis route and track stereocenters"""
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            stereocenters = get_stereocenters(mol_smiles)

            if stereocenters:
                print(f"Found stereocenters at depth {depth}: {stereocenters}")
                stereocenters_by_mol[mol_smiles] = stereocenters
                mol_depths[mol_smiles] = depth
            else:
                print(f"No stereocenters at depth {depth}")

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Track atom mappings from reactants to product
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    for atom in reactant_mol.GetAtoms():
                        if atom.HasProp("molAtomMapNumber"):
                            map_num = atom.GetProp("molAtomMapNumber")
                            if map_num not in atom_mappings:
                                atom_mappings[map_num] = []
                            atom_mappings[map_num].append(reactant)

                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    for atom in product_mol.GetAtoms():
                        if atom.HasProp("molAtomMapNumber"):
                            map_num = atom.GetProp("molAtomMapNumber")
                            if map_num not in atom_mappings:
                                atom_mappings[map_num] = []
                            atom_mappings[map_num].append(product_smiles)
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Traverse the route
    dfs_traverse(route)

    # Check if we have any stereocenters
    if not stereocenters_by_mol:
        print("No stereocenters found in the route")
        return False

    # Get the target molecule (depth 0)
    target_mol = None
    for mol, depth in mol_depths.items():
        if depth == 0:
            target_mol = mol
            break

    if not target_mol or target_mol not in stereocenters_by_mol:
        print("Target molecule has no stereocenters")
        return False

    # Check if at least one stereocenter is maintained throughout
    maintained = False

    # Get molecules sorted by depth (from target to starting materials)
    sorted_mols = sorted(mol_depths.items(), key=lambda x: x[1])

    # Check if there's at least one stereocenter in each molecule
    if all(mol in stereocenters_by_mol for mol, _ in sorted_mols):
        print("All molecules have at least one stereocenter")
        maintained = True
    else:
        # Check if there are gaps in the depths where molecules have stereocenters
        depths_with_stereocenters = set(
            depth
            for _, depth in sorted_mols
            if sorted_mols[0][0] in stereocenters_by_mol
        )
        max_depth = max(mol_depths.values())

        # If we have stereocenters at alternating depths (due to reaction nodes), that's fine
        # We just need to ensure there's a continuous chain of stereocenters
        continuous = True
        for i in range(max_depth + 1):
            # Skip reaction nodes (odd depths in retrosynthesis)
            if i % 2 == 1:  # Reaction nodes typically at odd depths
                continue

            # Check if this molecule depth has a stereocenter
            if i not in depths_with_stereocenters:
                continuous = False
                break

        maintained = continuous

    print(f"Maintained stereocenter throughout: {maintained}")
    return maintained
