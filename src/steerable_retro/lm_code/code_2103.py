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
    Detects if the synthetic route preserves a stereocenter throughout the synthesis.
    """
    # Track stereocenters by molecule SMILES and atom indices
    stereocenters_by_mol = {}
    # Track depth of each molecule
    depth_by_mol = {}
    # Track atom mappings through reactions
    atom_mappings = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            depth_by_mol[mol_smiles] = depth

            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                # Find chiral centers
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                if chiral_centers:
                    # Store atom indices of stereocenters
                    stereocenters_by_mol[mol_smiles] = set(
                        atom_idx for atom_idx, _ in chiral_centers
                    )
                    print(
                        f"Found {len(chiral_centers)} stereocenters at depth {depth} in molecule {mol_smiles}"
                    )

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract atom mappings from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            try:
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Create molecules from reactants and product
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Extract atom mapping numbers from product
                    product_atom_maps = {}
                    for atom in product_mol.GetAtoms():
                        map_num = atom.GetAtomMapNum()
                        if map_num > 0:
                            product_atom_maps[map_num] = atom.GetIdx()

                    # Store atom mappings for this reaction
                    atom_mappings[rsmi] = product_atom_maps
                    print(f"Extracted atom mappings from reaction: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Traverse the route to collect stereocenters and depths
    dfs_traverse(route)

    # Check if any stereocenter is preserved throughout the synthesis
    if len(stereocenters_by_mol) >= 2:
        # Sort molecules by depth
        mols_by_depth = sorted([(depth, smiles) for smiles, depth in depth_by_mol.items()])

        # Check if there's a significant depth difference (at least 2 levels)
        if len(mols_by_depth) >= 2:
            min_depth = mols_by_depth[0][0]
            max_depth = mols_by_depth[-1][0]

            if max_depth - min_depth >= 2:
                # Check if any stereocenter is preserved from early to late stages
                early_mols = [smiles for depth, smiles in mols_by_depth if depth >= max_depth - 1]
                late_mols = [smiles for depth, smiles in mols_by_depth if depth <= min_depth + 1]

                for early_mol in early_mols:
                    if early_mol in stereocenters_by_mol:
                        for late_mol in late_mols:
                            if late_mol in stereocenters_by_mol:
                                print(
                                    f"Stereocenter preservation detected from depth {max_depth} to {min_depth}"
                                )
                                print(f"Early molecule: {early_mol}")
                                print(f"Late molecule: {late_mol}")
                                return True

    return False
