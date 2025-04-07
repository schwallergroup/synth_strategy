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
    This function detects if the synthetic route contains both chloroaryl and benzodioxane motifs.
    """
    has_chloroaryl = False
    has_benzodioxane = False

    def dfs_traverse(node):
        nonlocal has_chloroaryl, has_benzodioxane

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for chloroaryl using the checker function
            if checker.check_fg("Aromatic halide", mol_smiles):
                # Verify it's specifically chlorine
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "Cl":
                            # Check if the chlorine is attached to an aromatic atom
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetIsAromatic():
                                    has_chloroaryl = True
                                    print(f"Found chloroaryl in molecule: {mol_smiles}")
                                    break
                            if has_chloroaryl:
                                break

            # Check for benzodioxane or related structures
            # Benzodioxane and its variants include structures where benzene is fused with
            # an oxygen-containing heterocycle like dioxane, dioxolane, or similar structures
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                # Check for benzene ring
                if checker.check_ring("benzene", mol_smiles):
                    # Check for various oxygen-containing rings that could form benzodioxane-like structures
                    oxygen_rings = ["dioxane", "dioxolane", "dioxolene"]
                    for ring in oxygen_rings:
                        if checker.check_ring(ring, mol_smiles):
                            has_benzodioxane = True
                            print(
                                f"Found benzodioxane-like structure with {ring} in molecule: {mol_smiles}"
                            )
                            break

                    # Also check for dihydrobenzodioxine pattern (benzene fused with a ring containing O and CH2 groups)
                    if not has_benzodioxane and "CCO" in mol_smiles:
                        # Look for patterns like c1ccc2c(c1)CCO2 which is dihydrobenzodioxine
                        dihydrobenzodioxine_pattern = Chem.MolFromSmarts("c1ccc2c(c1)CCO2")
                        if mol.HasSubstructMatch(dihydrobenzodioxine_pattern):
                            has_benzodioxane = True
                            print(f"Found dihydrobenzodioxine in molecule: {mol_smiles}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Final result - Found chloroaryl: {has_chloroaryl}, Found benzodioxane: {has_benzodioxane}"
    )

    return has_chloroaryl and has_benzodioxane
