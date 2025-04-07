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
    This function detects if the synthesis includes an aldol-type disconnection
    (breaking a β-hydroxy ketone into a ketone and an aldehyde).
    """
    has_aldol_disconnection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aldol_disconnection

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is an aldol condensation reaction
                if checker.check_reaction("Aldol condensation", rsmi):
                    print(f"Aldol condensation reaction detected at depth {depth}")
                    has_aldol_disconnection = True
                else:
                    # Alternative check: look for β-hydroxy ketone in product and ketone/aldehyde in reactants
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check if product has a β-hydroxy ketone structure
                    has_beta_hydroxy_ketone = False
                    if checker.check_fg("Ketone", product_smiles):
                        # Check for hydroxyl group in correct position relative to ketone
                        if product_mol:
                            beta_hydroxy_pattern = Chem.MolFromSmarts(
                                "[#6][C;$(C=O)][#6][C;!$(C=O)][OH]"
                            )
                            if product_mol.HasSubstructMatch(beta_hydroxy_pattern):
                                has_beta_hydroxy_ketone = True
                                print(f"β-hydroxy ketone found in product at depth {depth}")

                    # Check if reactants contain ketone and aldehyde
                    has_ketone = any(
                        checker.check_fg("Ketone", reactant) for reactant in reactants_smiles
                    )
                    has_aldehyde = any(
                        checker.check_fg("Aldehyde", reactant) for reactant in reactants_smiles
                    )

                    if has_beta_hydroxy_ketone and (has_ketone or has_aldehyde):
                        print(f"Aldol-type disconnection detected at depth {depth}")
                        print(f"Product: {product_smiles}")
                        print(f"Reactants: {reactants_smiles}")
                        has_aldol_disconnection = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_aldol_disconnection
