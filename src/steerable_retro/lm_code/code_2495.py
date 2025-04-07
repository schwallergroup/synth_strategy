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
    Detects if the synthesis route involves a transformation of a hydroxyl group
    to a chloro group, particularly on a heterocycle.
    """
    hydroxyl_to_chloro = False

    def dfs_traverse(node, depth=0):
        nonlocal hydroxyl_to_chloro

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                print(f"No rsmi found at depth {depth}")
                return

            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is an alcohol to chloride reaction
            alcohol_to_chloride_reactions = [
                "Alcohol to chloride_sulfonyl chloride",
                "Alcohol to chloride_SOCl2",
                "Alcohol to chloride_CHCl3",
                "Alcohol to chloride_CH2Cl2",
                "Alcohol to chloride_PCl5_ortho",
                "Alcohol to chloride_POCl3_ortho",
                "Alcohol to chloride_POCl3_para",
                "Alcohol to chloride_POCl3",
                "Alcohol to chloride_HCl",
                "Alcohol to chloride_Salt",
                "Alcohol to chloride_Other",
            ]

            is_alcohol_to_chloride = False
            for rxn_type in alcohol_to_chloride_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Detected {rxn_type} reaction at depth {depth}")
                    is_alcohol_to_chloride = True
                    break

            # Additional check for POCl3 or similar reagents in the reaction
            if not is_alcohol_to_chloride and (
                "O=P(Cl)(Cl)Cl" in rsmi or "ClP(Cl)(Cl)=O" in rsmi or "O=P(Cl)(Cl)[Cl" in rsmi
            ):
                print(f"Detected POCl3 reagent in reaction at depth {depth}")
                is_alcohol_to_chloride = True

            if is_alcohol_to_chloride or "O=P(Cl)(Cl)" in rsmi:
                # Check for heterocycles in both reactants and product
                heterocycles = [
                    "pyridine",
                    "pyrimidine",
                    "pyrazine",
                    "pyridazine",
                    "triazole",
                    "tetrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "furan",
                    "thiophene",
                    "pyrrole",
                    "isoxazole",
                    "isothiazole",
                ]

                # Check if product contains a heterocycle and a chloride
                product_has_heterocycle = any(
                    checker.check_ring(ring, product) for ring in heterocycles
                )
                product_has_chloride = "Cl" in product

                print(
                    f"Product has heterocycle: {product_has_heterocycle}, Product has chloride: {product_has_chloride}"
                )

                if product_has_heterocycle and product_has_chloride:
                    # Check if any reactant has a hydroxyl group and a heterocycle
                    for reactant in reactants:
                        reactant_has_heterocycle = any(
                            checker.check_ring(ring, reactant) for ring in heterocycles
                        )

                        # Check for any type of alcohol
                        alcohol_types = [
                            "Phenol",
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol",
                            "Aromatic alcohol",
                        ]
                        reactant_has_alcohol = any(
                            checker.check_fg(alcohol, reactant) for alcohol in alcohol_types
                        )

                        # Also check for OH group directly
                        reactant_has_oh = (
                            "OH" in reactant or "O[H]" in reactant or "[OH]" in reactant
                        )

                        print(f"Reactant: {reactant}")
                        print(f"  Has heterocycle: {reactant_has_heterocycle}")
                        print(f"  Has alcohol: {reactant_has_alcohol}")
                        print(f"  Has OH group: {reactant_has_oh}")

                        if reactant_has_heterocycle and (reactant_has_alcohol or reactant_has_oh):
                            # We've found a reaction that converts a hydroxyl on a heterocycle to a chloro
                            print(
                                f"Found hydroxyl to chloro transformation on heterocycle at depth {depth}"
                            )
                            hydroxyl_to_chloro = True
                            return  # Found what we're looking for, no need to continue

            # Special case for heterocyclic OH to Cl transformation
            if "Oc" in reactants_part and "Cl" in product and not is_alcohol_to_chloride:
                for heterocycle in [
                    "pyridine",
                    "pyrimidine",
                    "pyrazine",
                    "pyridazine",
                    "triazole",
                    "tetrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                ]:
                    for reactant in reactants:
                        if checker.check_ring(heterocycle, reactant) and (
                            "OH" in reactant
                            or "O[H]" in reactant
                            or "[OH]" in reactant
                            or "Oc" in reactant
                        ):
                            if checker.check_ring(heterocycle, product) and "Cl" in product:
                                print(
                                    f"Found direct heterocyclic OH to Cl transformation at depth {depth}"
                                )
                                hydroxyl_to_chloro = True
                                return

        # Continue traversing
        for child in node.get("children", []):
            if not hydroxyl_to_chloro:  # Only continue if we haven't found the transformation yet
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {hydroxyl_to_chloro}")
    return hydroxyl_to_chloro
