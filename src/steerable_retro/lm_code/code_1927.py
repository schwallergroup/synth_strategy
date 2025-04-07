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
    This function detects the use of Williamson ether synthesis in the route,
    specifically looking for the reaction between a benzyl halide and a phenol.
    """
    has_williamson_ether = False

    def dfs_traverse(node):
        nonlocal has_williamson_ether

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # First check if this is a Williamson ether synthesis reaction by name
                if checker.check_reaction(
                    "Williamson Ether Synthesis", rsmi
                ) or checker.check_reaction("{Williamson ether}", rsmi):
                    print(f"Detected Williamson ether synthesis by reaction name: {rsmi}")
                    has_williamson_ether = True

                # Alternative check for specific case of benzyl halide + phenol
                elif not has_williamson_ether:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reaction: {rsmi}")

                    # Check for phenol in reactants - expanded definition
                    has_phenol = any(
                        checker.check_fg("Phenol", reactant)
                        or (
                            ("[OH]" in reactant or "OH" in reactant)
                            and checker.check_ring("benzene", reactant)
                        )
                        for reactant in reactants
                    )
                    print(f"Has phenol: {has_phenol}")

                    # Check for benzyl halide (halogen attached to CH2 connected to benzene)
                    has_benzyl_halide = False
                    benzyl_halide_reactant = None
                    for reactant in reactants:
                        # Check for primary halide attached to benzene (benzyl halide)
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        ) and checker.check_ring("benzene", reactant):
                            # Additional check for CH2 group
                            if "[CH2]" in reactant or "CH2" in reactant:
                                has_benzyl_halide = True
                                benzyl_halide_reactant = reactant
                                print(f"Found benzyl halide: {reactant}")
                                break

                    # Check if product has ether and benzene ring
                    has_ether = checker.check_fg("Ether", product)
                    has_benzene_in_product = checker.check_ring("benzene", product)
                    print(
                        f"Has ether in product: {has_ether}, Has benzene in product: {has_benzene_in_product}"
                    )

                    # Check for the specific reaction pattern
                    if has_phenol and has_benzyl_halide and has_ether and has_benzene_in_product:
                        print(
                            f"Detected Williamson ether synthesis (benzyl halide + phenol): {rsmi}"
                        )
                        has_williamson_ether = True
                    # Also check for general Williamson pattern (any halide + phenol)
                    elif (
                        has_phenol
                        and any(
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            for reactant in reactants
                        )
                        and has_ether
                    ):
                        print(f"Detected general Williamson ether synthesis: {rsmi}")
                        has_williamson_ether = True
                    # Check for benzyl halide + any OH group forming an ether
                    elif (
                        has_benzyl_halide
                        and has_ether
                        and has_benzene_in_product
                        and any("[OH]" in reactant or "OH" in reactant for reactant in reactants)
                    ):
                        print(
                            f"Detected Williamson-like ether synthesis (benzyl halide + OH): {rsmi}"
                        )
                        has_williamson_ether = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_williamson_ether
