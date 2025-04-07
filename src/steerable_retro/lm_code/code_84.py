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
    This function detects if the synthetic route involves sequential functional group transformations
    (nitration followed by reduction).
    """
    # Track reactions by depth
    nitration_depth = None
    reduction_depth = None

    def dfs_traverse(node, current_depth=0):
        nonlocal nitration_depth, reduction_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitration reaction
                if (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                ):
                    print(f"Nitration detected at depth {current_depth}")
                    nitration_depth = current_depth

                # Check for reduction of nitro to amine
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Reduction detected at depth {current_depth}")
                    reduction_depth = current_depth

                # If specific reaction types aren't detected, check for functional group changes
                if nitration_depth is None:
                    # Check if any reactant doesn't have nitro group but product does
                    product_has_nitro = checker.check_fg("Nitro group", product_smiles)
                    reactants_all_have_nitro = all(
                        checker.check_fg("Nitro group", r) for r in reactants_smiles
                    )

                    if product_has_nitro and not reactants_all_have_nitro:
                        print(f"Nitration (FG change) detected at depth {current_depth}")
                        nitration_depth = current_depth

                if reduction_depth is None:
                    # Check if any reactant has nitro group but product has amine instead
                    reactant_has_nitro = any(
                        checker.check_fg("Nitro group", r) for r in reactants_smiles
                    )
                    product_has_amine = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )

                    if reactant_has_nitro and product_has_amine:
                        print(f"Reduction (FG change) detected at depth {current_depth}")
                        reduction_depth = current_depth

        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    # Check if nitration is followed by reduction
    if nitration_depth is not None and reduction_depth is not None:
        # In retrosynthetic direction, higher depth is earlier in synthesis
        if nitration_depth > reduction_depth:
            print("Sequential nitration-reduction strategy detected")
            return True

    return False
