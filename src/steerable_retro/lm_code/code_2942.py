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

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects if aryl chlorides are maintained throughout the synthesis.
    This checks if aryl chlorides present in the target molecule are preserved
    in all steps of the synthesis route.
    """
    # Track if we've found an aryl chloride in the target molecule
    target_has_aryl_chloride = False
    # Track if any step loses the aryl chloride
    aryl_chloride_maintained = True

    def has_aromatic_chloride(mol_smiles):
        """Helper function to check if a molecule has an aromatic chloride"""
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol is None:
                return False

            # Check if molecule has aromatic halide functional group
            if not checker.check_fg("Aromatic halide", mol_smiles):
                return False

            # Verify that at least one of the aromatic halides is a chlorine
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "Cl" and atom.GetNeighbors()[0].GetIsAromatic():
                    return True
            return False
        except:
            # Fallback to simpler check if RDKit parsing fails
            return checker.check_fg("Aromatic halide", mol_smiles) and "Cl" in mol_smiles

    def dfs_traverse(node, depth=0, parent_has_aryl_cl=None):
        nonlocal target_has_aryl_chloride, aryl_chloride_maintained

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_aryl_chloride = has_aromatic_chloride(mol_smiles)

            # If this is the target molecule (depth 0)
            if depth == 0:
                target_has_aryl_chloride = has_aryl_chloride
                print(f"Target molecule has aryl chloride: {target_has_aryl_chloride}")

                # If target doesn't have aryl chloride, no need to check maintenance
                if not target_has_aryl_chloride:
                    return

            # For non-target molecules, check if aryl chloride is maintained
            elif parent_has_aryl_cl is not None:
                # Only check non-starting materials
                if parent_has_aryl_cl and not has_aryl_chloride and not node.get("in_stock", False):
                    print(f"Aryl chloride lost at depth {depth} in molecule: {mol_smiles}")
                    aryl_chloride_maintained = False

            # Pass down whether this molecule has aryl chloride
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, has_aryl_chloride)

        elif node["type"] == "reaction":
            # For reaction nodes, we just pass through the parent's aryl chloride status
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, parent_has_aryl_cl)

    # Start traversal
    dfs_traverse(route)

    # If target doesn't have aryl chloride, return False (no maintenance to check)
    if not target_has_aryl_chloride:
        print("Target molecule doesn't have aryl chloride, so maintenance check is not applicable")
        return False

    print(f"Aryl chloride maintained throughout synthesis: {aryl_chloride_maintained}")
    return aryl_chloride_maintained
