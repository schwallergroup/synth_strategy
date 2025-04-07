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
    This function detects if the synthesis involves ketone protection as a ketal.
    """
    protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # In retrosynthesis, we're looking for ketal hydrolysis (deprotection)
            # which would be the reverse of protection in forward synthesis

            # Check if product has ketone
            product_has_ketone = checker.check_fg("Ketone", product)

            # Check if any reactant contains ketal
            reactant_has_ketal = any(
                checker.check_fg("Acetal/Ketal", reactant) for reactant in reactants
            )

            # Check if this is a ketal hydrolysis reaction
            if checker.check_reaction(
                "Acetal hydrolysis to ketone", rsmi
            ) or checker.check_reaction("Ketal hydrolysis to ketone", rsmi):
                print(f"Detected ketal deprotection to ketone: {rsmi}")
                protection_found = True
                return

            # Alternative check: reactant has ketal, product has ketone, and no ketal in product
            product_has_ketal = checker.check_fg("Acetal/Ketal", product)
            if reactant_has_ketal and product_has_ketone and not product_has_ketal:
                print(f"Detected ketal deprotection to ketone: {rsmi}")
                protection_found = True
                return

        # Continue DFS traversal
        for child in node.get("children", []):
            if not protection_found:  # Stop traversal if protection already found
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return protection_found
