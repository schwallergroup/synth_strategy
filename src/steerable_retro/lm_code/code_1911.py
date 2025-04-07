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
    Detects if the synthesis involves an azo coupling strategy where a diazonium
    intermediate couples with an aromatic compound.
    """
    has_azo_formation = False

    def dfs_traverse(node):
        nonlocal has_azo_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains azo group (diazene)
            if checker.check_fg("Diazene", product):
                print(f"Found product with diazene (azo) group: {product}")

                # Check if reactants contain compounds that can participate in azo coupling
                for reactant in reactants:
                    # Check for diazonium precursors or intermediates
                    if (
                        checker.check_fg("Diazo", reactant)
                        or checker.check_fg("Azide", reactant)
                        or checker.check_fg("Nitroso", reactant)
                        or checker.check_fg("Nitro group", reactant)
                        or checker.check_fg("Aniline", reactant)
                        or "N=[N+]" in reactant
                    ):  # Explicit check for diazonium salt

                        print(f"Found reactant that can participate in azo coupling: {reactant}")
                        has_azo_formation = True
                        break

            # Check if this is a diazotization reaction (forming diazonium salt)
            if not has_azo_formation and "N=[N+]" in product:
                print(f"Found potential diazonium formation: {rsmi}")
                for reactant in reactants:
                    if checker.check_fg("Aniline", reactant) and (
                        checker.check_fg("Nitroso", reactant)
                        or checker.check_fg("Nitro group", reactant)
                        or "NO2" in reactant
                        or "NO3" in reactant
                        or "ONO" in reactant
                    ):
                        print(f"Found diazotization reaction with aniline: {reactant}")
                        has_azo_formation = True
                        break

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_azo_formation
