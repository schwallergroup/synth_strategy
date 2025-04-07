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
    Detects if the synthetic route preserves stereocenters throughout the synthesis.

    Returns True if:
    1. The final product has at least one stereocenter
    2. All reactions in the route preserve existing stereocenters
    """
    # Track if we've found any violations of stereocenter preservation
    stereocenters_preserved = True
    # Track if the final product has stereocenters
    has_final_product_stereocenters = False

    def get_chiral_atoms_with_mapping(mol_smiles):
        """Helper function to get chiral atoms with their atom mapping"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return {}

        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
        chiral_atoms = {}

        for atom_idx, chirality in chiral_centers:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Get atom mapping if available
            atom_map = (
                atom.GetProp("molAtomMapNumber") if atom.HasProp("molAtomMapNumber") else None
            )
            if atom_map:
                chiral_atoms[atom_map] = chirality

        return chiral_atoms

    def dfs_traverse(node, depth=0):
        nonlocal stereocenters_preserved, has_final_product_stereocenters

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and depth == 0:  # Final product
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    has_final_product_stereocenters = len(chiral_centers) > 0
                    print(f"Final product has {len(chiral_centers)} stereocenters")
            except Exception as e:
                print(f"Error processing molecule: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Get chiral atoms with their mapping from reactants and product
                reactant_chiral_atoms = get_chiral_atoms_with_mapping(reactants_smiles)
                product_chiral_atoms = get_chiral_atoms_with_mapping(product_smiles)

                # Check if any chiral atoms in reactants lost their chirality in the product
                for atom_map, chirality in reactant_chiral_atoms.items():
                    if atom_map in product_chiral_atoms:
                        # Check if chirality is preserved (R->R or S->S)
                        if product_chiral_atoms[atom_map] != chirality:
                            print(
                                f"Stereocenter with map {atom_map} changed from {chirality} to {product_chiral_atoms[atom_map]}"
                            )
                            stereocenters_preserved = False
                    else:
                        # Mapped atom exists in reactants but not in products or lost chirality
                        print(f"Stereocenter with map {atom_map} was lost in the product")
                        stereocenters_preserved = False
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return stereocenters_preserved and has_final_product_stereocenters
