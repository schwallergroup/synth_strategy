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
    This function detects preservation of cyano groups throughout the synthesis.
    In retrosynthesis, we check if cyano groups in the product are preserved in the reactants.
    """
    cyano_groups_preserved = True
    reactions_with_cyano = 0

    def dfs_traverse(node, depth=0):
        nonlocal cyano_groups_preserved, reactions_with_cyano

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has cyano groups
                product_has_cyano = checker.check_fg("Nitrile", product_smiles)

                if product_has_cyano:
                    reactions_with_cyano += 1

                    # In retrosynthesis, we check if all cyano groups in product are preserved in reactants
                    reactants_have_all_cyano = False

                    # Get atom indices of cyano groups in product
                    product_cyano_indices = checker.get_fg_atom_indices("Nitrile", product_smiles)

                    if product_cyano_indices:
                        # Check if all product cyano groups are preserved in reactants
                        # This is a simplification - ideally we would use atom mapping to track specific groups
                        reactant_has_cyano = False
                        for reactant in reactants_smiles:
                            if checker.check_fg("Nitrile", reactant):
                                reactant_has_cyano = True
                                break

                        if not reactant_has_cyano:
                            cyano_groups_preserved = False
                            print(
                                f"Cyano group not preserved at depth {depth}: product has cyano but reactants don't"
                            )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Only return True if we actually had cyano groups and they were preserved
    return cyano_groups_preserved and reactions_with_cyano > 0
