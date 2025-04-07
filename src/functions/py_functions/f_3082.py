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
    This function detects if the synthetic route involves nitration of an aromatic amine.
    """
    has_nitration = False

    def dfs_traverse(node):
        nonlocal has_nitration

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                try:
                    # Correctly split reaction SMILES
                    reactants = rsmi.split(">")[0]
                    product = rsmi.split(">")[-1]

                    # Check for aromatic nitration reaction
                    is_nitration_reaction = (
                        checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                        or checker.check_reaction(
                            "Aromatic nitration with NO3 salt", rsmi
                        )
                        or checker.check_reaction(
                            "Aromatic nitration with NO2 salt", rsmi
                        )
                        or checker.check_reaction(
                            "Aromatic nitration with alkyl NO2", rsmi
                        )
                    )

                    # Check for aromatic ring in reactants
                    for reactant in reactants.split("."):
                        # Check if reactant has an aromatic ring and an amine group
                        has_aromatic = (
                            checker.check_ring("benzene", reactant)
                            or checker.check_ring("naphthalene", reactant)
                            or checker.check_ring("anthracene", reactant)
                            or checker.check_ring("pyridine", reactant)
                            or checker.check_ring("pyrrole", reactant)
                            or checker.check_ring("furan", reactant)
                            or checker.check_ring("thiophene", reactant)
                        )

                        has_amine = (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        )

                        # Check if product has a nitro group
                        has_nitro = checker.check_fg("Nitro group", product)

                        if (
                            is_nitration_reaction
                            and has_aromatic
                            and has_amine
                            and has_nitro
                        ):
                            print(f"Found nitration of aromatic amine: {rsmi}")
                            has_nitration = True
                            break

                        # Also check for direct transformation of amine to nitro group
                        # This would be a different reaction than aromatic nitration
                        if (
                            has_amine
                            and has_nitro
                            and not checker.check_fg("Aniline", product)
                        ):
                            print(
                                f"Found transformation of amine to nitro group: {rsmi}"
                            )
                            has_nitration = True
                            break
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Nitro group transformation detected: {has_nitration}")
    return has_nitration
