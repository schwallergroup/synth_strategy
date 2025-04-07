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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the synthesis route includes nitro reduction to amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction" and not nitro_reduction_detected:
            try:
                # Get reaction SMILES
                rsmi = node["metadata"]["rsmi"]

                # Direct check for nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print("Detected nitro reduction to amine via reaction pattern")
                    nitro_reduction_detected = True
                    return

                # Fallback: Check for nitro to amine conversion manually
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has nitro group and product has amine
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant):
                        # Found a reactant with nitro group, now check if product has amine
                        if (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        ):
                            # Additional check to ensure it's a reduction reaction
                            # Look for reducing conditions or hydrogen addition
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            product_mol = Chem.MolFromSmiles(product)

                            if reactant_mol and product_mol:
                                # Check if the molecule with nitro has fewer oxygens in the product
                                # which would indicate reduction
                                reactant_o_count = sum(
                                    1
                                    for atom in reactant_mol.GetAtoms()
                                    if atom.GetSymbol() == "O"
                                )
                                product_o_count = sum(
                                    1
                                    for atom in product_mol.GetAtoms()
                                    if atom.GetSymbol() == "O"
                                )

                                if product_o_count < reactant_o_count:
                                    print(
                                        f"Detected nitro reduction to amine via functional group analysis"
                                    )
                                    print(f"Reactant: {reactant}")
                                    print(f"Product: {product}")
                                    nitro_reduction_detected = True
                                    return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return nitro_reduction_detected
