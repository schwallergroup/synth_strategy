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
    Detects if the synthesis introduces a nitrile-containing fragment through a key disconnection.
    """
    nitrile_introduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_introduction_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Skip if we already found a nitrile introduction
            if nitrile_introduction_found:
                return

            # Check if product has nitrile
            product_has_nitrile = checker.check_fg("Nitrile", product)

            if product_has_nitrile and len(reactants) >= 2:
                # Check if some reactants have nitrile but not all
                reactants_with_nitrile = [r for r in reactants if checker.check_fg("Nitrile", r)]

                if len(reactants_with_nitrile) > 0 and len(reactants_with_nitrile) < len(reactants):
                    # Check if this is a coupling reaction (common for fragment introduction)
                    is_coupling = any(
                        [
                            checker.check_reaction("Suzuki coupling with boronic acids", rsmi),
                            checker.check_reaction("Suzuki coupling with boronic esters", rsmi),
                            checker.check_reaction("Negishi coupling", rsmi),
                            checker.check_reaction("Stille reaction_aryl", rsmi),
                            checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi),
                        ]
                    )

                    if is_coupling:
                        print(
                            f"Found nitrile-containing fragment introduction through coupling at depth {depth}"
                        )
                        nitrile_introduction_found = True
                    else:
                        # Even if not a known coupling, it might still be a fragment introduction
                        print(
                            f"Found potential nitrile-containing fragment introduction at depth {depth}"
                        )
                        nitrile_introduction_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitrile_introduction_found
