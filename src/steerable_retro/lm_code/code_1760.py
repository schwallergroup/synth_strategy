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
    This function detects the transformation of an oxime to an oxime ether.
    This is typically done via a Williamson ether synthesis reaction.
    """
    oxime_ether_formation_found = False

    def dfs_traverse(node):
        nonlocal oxime_ether_formation_found

        if node["type"] == "reaction" and not oxime_ether_formation_found:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if this is a Williamson ether synthesis
                    is_williamson = checker.check_reaction("Williamson Ether Synthesis", rsmi)

                    # Check for oxime in reactants
                    has_oxime_reactant = any(checker.check_fg("Oxime", r) for r in reactants)

                    # Check for primary/secondary/tertiary halide in reactants
                    has_halide_reactant = any(
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        for r in reactants
                    )

                    # Check for oxime ether in product
                    # Since "Oxime ether" might not be in the FG list, we check for absence of oxime
                    # and presence of ether in the product
                    has_oxime_product = checker.check_fg("Oxime", product)
                    has_ether_product = checker.check_fg("Ether", product)

                    if (
                        is_williamson
                        and has_oxime_reactant
                        and has_halide_reactant
                        and not has_oxime_product
                        and has_ether_product
                    ):
                        print(
                            "Found oxime to oxime ether transformation via Williamson ether synthesis"
                        )
                        oxime_ether_formation_found = True

                    # Fallback check if Williamson reaction check fails
                    elif has_oxime_reactant and has_halide_reactant and not has_oxime_product:
                        # Use SMARTS pattern for oxime ether as fallback
                        oxime_ether_pattern = Chem.MolFromSmarts("[C]=[N][O][C]")
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(oxime_ether_pattern):
                            print("Found oxime to oxime ether transformation (fallback detection)")
                            oxime_ether_formation_found = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return oxime_ether_formation_found
