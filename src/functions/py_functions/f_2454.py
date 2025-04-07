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
    This function detects the use of phthalimide for amine protection in a synthesis route.
    It looks for a sequence where a primary amine is protected with phthalimide.
    """
    phthalimide_used = False

    def dfs_traverse(node):
        nonlocal phthalimide_used

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phthalimide protection reaction directly
                if checker.check_reaction("Phthalic anhydride to phthalimide", rsmi):
                    print(f"Found phthalimide formation reaction: {rsmi}")

                    # Verify this is for protection by checking if primary amine is involved
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant):
                            print(
                                f"Found primary amine in phthalimide formation: {reactant}"
                            )
                            phthalimide_used = True
                            return

                # Check for primary amine in reactants
                primary_amine_reactant = None
                for reactant in reactants:
                    if checker.check_fg("Primary amine", reactant):
                        primary_amine_reactant = reactant
                        print(f"Found primary amine in reactants: {reactant}")
                        break

                if primary_amine_reactant:
                    # Check if product contains a protected amine structure
                    # This would typically be a secondary or tertiary amide in a cyclic structure
                    if checker.check_fg("Secondary amide", product) or checker.check_fg(
                        "Tertiary amide", product
                    ):

                        # Check if the product has a phthalimide-like structure
                        if checker.check_ring("phthalimide", product):
                            print(f"Found phthalimide protection reaction: {rsmi}")
                            phthalimide_used = True
                            return

                        # Check for other cyclic imide structures that might be phthalimide-related
                        for ring_type in ["benzene", "naphthalene"]:
                            if checker.check_ring(ring_type, product):
                                print(
                                    f"Found potential phthalimide-like protection with {ring_type}: {rsmi}"
                                )
                                phthalimide_used = True
                                return

                # Check for phthalimide deprotection (should not count as protection)
                if checker.check_reaction("Phthalimide deprotection", rsmi):
                    print(
                        f"Found phthalimide deprotection (not counted as protection): {rsmi}"
                    )
                    return

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return phthalimide_used
