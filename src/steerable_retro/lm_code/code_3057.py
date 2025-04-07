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
    This function detects a synthetic strategy involving early-stage halogenation,
    specifically looking for halogenation reactions (bromination, chlorination,
    fluorination, iodination) at high depth values (early in synthesis).
    """
    early_halogenation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal early_halogenation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an early stage reaction (depth >= 3)
                if current_depth >= 3:
                    # Check for halogenation reactions using the checker function
                    halogenation_reactions = [
                        "Aromatic bromination",
                        "Aromatic chlorination",
                        "Aromatic fluorination",
                        "Aromatic iodination",
                        "Bromination",
                        "Chlorination",
                        "Fluorination",
                        "Iodination",
                    ]

                    for reaction_type in halogenation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Early-stage {reaction_type} detected at depth {current_depth}")
                            early_halogenation = True
                            break

                    # If no specific halogenation reaction was detected, check for halogen addition
                    if not early_halogenation:
                        product_mol = Chem.MolFromSmiles(product)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r.strip()]

                        # Check if product has halogen that wasn't in reactants
                        has_new_halogen = False
                        halogen_fgs = [
                            "Primary halide",
                            "Secondary halide",
                            "Tertiary halide",
                            "Aromatic halide",
                            "Alkenyl halide",
                            "Haloalkyne",
                        ]

                        for fg in halogen_fgs:
                            if checker.check_fg(fg, product):
                                # Check if this halogen was newly introduced
                                if not all(checker.check_fg(fg, r) for r in reactants if r.strip()):
                                    print(
                                        f"Early-stage halogenation (adding {fg}) detected at depth {current_depth}"
                                    )
                                    early_halogenation = True
                                    break
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    return early_halogenation
