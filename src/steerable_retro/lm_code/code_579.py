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
    This function detects a synthetic strategy involving formation of a heterocycle
    with two phenyl substituents.
    """
    has_heterocycle_formation = False
    has_diphenyl_substitution = False

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation, has_diphenyl_substitution

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocycle formation
            heterocycle_types = [
                "oxazole",
                "pyrazole",
                "imidazole",
                "thiazole",
                "pyrimidine",
                "isoxazole",
                "isothiazole",
                "triazole",
                "tetrazole",
                "furan",
                "pyrrole",
                "thiophene",
                "pyridine",
            ]

            # Check if product has a heterocycle
            product_heterocycle = None
            for htype in heterocycle_types:
                if checker.check_ring(htype, product_smiles):
                    product_heterocycle = htype
                    break

            # Check if any reactant has the same heterocycle
            reactant_has_same_heterocycle = False
            if product_heterocycle:
                for r_smiles in reactants_smiles:
                    if checker.check_ring(product_heterocycle, r_smiles):
                        reactant_has_same_heterocycle = True
                        break

            # Heterocycle formation occurs when product has heterocycle but reactants don't
            heterocycle_formed = product_heterocycle and not reactant_has_same_heterocycle

            # Check for benzene rings in product
            benzene_count = 0
            if checker.check_ring("benzene", product_smiles):
                # Get all benzene ring indices
                benzene_indices = checker.get_ring_atom_indices("benzene", product_smiles)
                benzene_count = len(benzene_indices)

            # Both conditions met in the same reaction
            if heterocycle_formed and benzene_count >= 2:
                has_heterocycle_formation = True
                has_diphenyl_substitution = True
                print(f"Found heterocycle formation with diphenyl substitution at depth {depth}")
                print(
                    f"Heterocycle type: {product_heterocycle}, Number of benzene rings: {benzene_count}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = has_heterocycle_formation and has_diphenyl_substitution

    print(f"Heterocycle formation with diphenyl substitution strategy detected: {strategy_present}")
    return strategy_present
