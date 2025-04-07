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
    This function detects a synthetic strategy involving late-stage formation of a tertiary alcohol
    via addition of a heterocyclic reagent to a ketone.
    """
    # Track if we found the tertiary alcohol formation
    found_tertiary_alcohol_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_tertiary_alcohol_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            product = product_part

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for tertiary alcohol in product
                has_tertiary_alcohol_in_product = checker.check_fg("Tertiary alcohol", product)
                print(f"Has tertiary alcohol in product: {has_tertiary_alcohol_in_product}")

                # Check for ketone in reactants
                has_ketone_in_reactants = False
                ketone_reactant = None
                for reactant in reactants:
                    if checker.check_fg("Ketone", reactant):
                        has_ketone_in_reactants = True
                        ketone_reactant = reactant
                        print(f"Found ketone in reactant: {reactant}")
                        break

                # Check for heterocyclic structure in reactants
                has_heterocycle = False
                heterocycle_reactant = None
                heterocycle_rings = [
                    "pyridine",
                    "furan",
                    "thiophene",
                    "pyrrole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "pyrimidine",
                    "pyrazine",
                    "pyridazine",
                    "triazole",
                    "tetrazole",
                    "indole",
                    "benzimidazole",
                    "quinoline",
                    "isoquinoline",
                    "benzoxazole",
                    "benzothiazole",
                ]

                for reactant in reactants:
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, reactant):
                            has_heterocycle = True
                            heterocycle_reactant = reactant
                            print(f"Found heterocycle ({ring}) in reactant: {reactant}")
                            break
                    if has_heterocycle:
                        break

                # Check if this is a Grignard or organolithium addition
                is_grignard_reaction = checker.check_reaction(
                    "Grignard from ketone to alcohol", rsmi
                )

                # Check for organolithium reaction - look for Li in any reactant
                is_organolithium_reaction = False
                for reactant in reactants:
                    if "Li" in reactant and has_ketone_in_reactants:
                        is_organolithium_reaction = True
                        print(f"Found organolithium reagent: {reactant}")
                        break

                # Check for general addition to ketone forming tertiary alcohol
                is_addition_to_ketone = False
                if has_tertiary_alcohol_in_product and has_ketone_in_reactants:
                    # Check that tertiary alcohol wasn't already in reactants
                    tertiary_alcohol_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Tertiary alcohol", reactant):
                            tertiary_alcohol_in_reactants = True
                            break

                    if not tertiary_alcohol_in_reactants:
                        is_addition_to_ketone = True
                        print("Confirmed tertiary alcohol formation from ketone")

                # If all conditions are met, we found the pattern
                if (
                    is_addition_to_ketone
                    and has_heterocycle
                    and (is_grignard_reaction or is_organolithium_reaction or "Li" in rsmi)
                ):
                    print(f"Found late-stage tertiary alcohol formation via addition to ketone")
                    print(f"Ketone reactant: {ketone_reactant}")
                    print(f"Heterocyclic reactant: {heterocycle_reactant}")
                    print(f"Product with tertiary alcohol: {product}")
                    found_tertiary_alcohol_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_tertiary_alcohol_formation
