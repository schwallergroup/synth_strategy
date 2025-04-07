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
    This function detects N-methylation of amines in the synthetic route.
    """
    has_n_methylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_n_methylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for methylation reactions directly
                methylation_reactions = [
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                    "Parnes methylation",
                    "N-methylation",
                    "Methylation with MeI_primary",
                    "Methylation with MeI_secondary",
                    "Methylation with MeI_tertiary",
                    "Methylation with DMS",
                    "DMS Amine methylation",
                    "Reductive amination with aldehyde",
                ]

                for reaction_type in methylation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found {reaction_type} at depth {depth}")
                        has_n_methylation = True
                        return

                # If no specific reaction type matched, check for amine methylation manually
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has an amine
                if (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                ):

                    # Check for formaldehyde in reactants
                    for reactant in reactants:
                        if checker.check_fg("Formaldehyde", reactant):
                            print(f"Found formaldehyde in reactants at depth {depth}")

                            # Check for amine conversion
                            for r in reactants:
                                if checker.check_fg("Primary amine", r) and checker.check_fg(
                                    "Secondary amine", product
                                ):
                                    print(
                                        f"Found primary to secondary amine methylation at depth {depth}"
                                    )
                                    has_n_methylation = True
                                    return
                                elif checker.check_fg("Secondary amine", r) and checker.check_fg(
                                    "Tertiary amine", product
                                ):
                                    print(
                                        f"Found secondary to tertiary amine methylation at depth {depth}"
                                    )
                                    has_n_methylation = True
                                    return

                    # Use atom mapping to detect methylation
                    try:
                        # Look for a methyl group (CH3) in the reactants that gets attached to an N in the product
                        for reactant in reactants:
                            # Check for formaldehyde or other C1 fragments
                            if (
                                "[CH2:15]" in reactant
                                or "[CH3:15]" in reactant
                                or "O=[CH2:15]" in reactant
                            ):
                                # Check if this carbon is attached to nitrogen in product
                                if "[N:14]([CH3:15])" in product or "[NH:14]([CH3:15])" in product:
                                    print(f"Found methylation via atom mapping at depth {depth}")
                                    has_n_methylation = True
                                    return
                    except Exception as e:
                        print(f"Error in atom mapping analysis: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_n_methylation
