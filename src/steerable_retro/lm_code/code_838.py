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
    Detects if the synthesis uses a late-stage cyanation strategy with convergent
    fragment assembly via sequential amide couplings.
    """
    has_late_cyanation = False
    has_amide_couplings = 0
    has_nitro_reduction = False
    has_acid_activation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_cyanation, has_amide_couplings, has_nitro_reduction, has_acid_activation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reaction_smiles = (
                    node["metadata"]["smiles"] if "smiles" in node["metadata"] else rsmi
                )
                product_smiles = rsmi.split(">")[-1]
                reactants_smiles = rsmi.split(">")[0].split(".")

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for cyanation at depth 0 or 1 (late-stage)
                if depth <= 1:
                    # Check for aromatic halide to nitrile conversion
                    if checker.check_fg("Nitrile", product_smiles) and any(
                        checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                    ):
                        print(f"Found late-stage cyanation at depth {depth}")
                        has_late_cyanation = True

                # Check for amide formation reactions - named reactions
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        reaction_smiles,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        reaction_smiles,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with secondary amine to amide", reaction_smiles
                    )
                    or checker.check_reaction(
                        "Carboxylic acid with primary amine to amide", reaction_smiles
                    )
                    or checker.check_reaction("Ester with primary amine to amide", reaction_smiles)
                    or checker.check_reaction(
                        "Ester with secondary amine to amide", reaction_smiles
                    )
                    or checker.check_reaction("Acylation of primary amines", reaction_smiles)
                    or checker.check_reaction("Acylation of secondary amines", reaction_smiles)
                    or checker.check_reaction("Schotten-Baumann to ester", reaction_smiles)
                    or checker.check_reaction("{Schotten-Baumann_amide}", reaction_smiles)
                ):
                    print(f"Found amide coupling (named reaction) at depth {depth}")
                    has_amide_couplings += 1

                # Check for amide formation by examining product and reactants
                elif (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                ):
                    # Check if reactants have amine and carboxylic acid/acyl halide
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants_smiles
                    )
                    has_acid = any(
                        checker.check_fg("Carboxylic acid", r)
                        or checker.check_fg("Acyl halide", r)
                        or checker.check_fg("Ester", r)
                        for r in reactants_smiles
                    )

                    if has_amine and has_acid:
                        print(f"Found amide coupling (FG analysis) at depth {depth}")
                        has_amide_couplings += 1

                # Check for nitro reduction
                if checker.check_reaction("Reduction of nitro groups to amines", reaction_smiles):
                    print(f"Found nitro reduction at depth {depth}")
                    has_nitro_reduction = True
                elif checker.check_fg(
                    "Nitro group", "".join(reactants_smiles)
                ) and checker.check_fg("Primary amine", product_smiles):
                    print(f"Found nitro reduction (FG analysis) at depth {depth}")
                    has_nitro_reduction = True

                # Check for acid activation
                if (
                    checker.check_reaction("Acyl chlorides from alcohols", reaction_smiles)
                    or checker.check_reaction("Acyl bromides from alcohols", reaction_smiles)
                    or checker.check_reaction("Acyl iodides from alcohols", reaction_smiles)
                ):
                    print(f"Found acid activation at depth {depth}")
                    has_acid_activation = True

                # Additional check for acid chloride formation from carboxylic acid
                if checker.check_fg("Acyl halide", product_smiles):
                    for reactant in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", reactant):
                            print(f"Found acid activation (COOH → COX) at depth {depth}")
                            has_acid_activation = True
                            break

                # Check for other acid activation methods
                for reactant in reactants_smiles:
                    if checker.check_fg("Carboxylic acid", reactant):
                        # Check if product has an activated form of the acid
                        if checker.check_fg("Ester", product_smiles) or checker.check_fg(
                            "Anhydride", product_smiles
                        ):
                            print(f"Found acid activation (other method) at depth {depth}")
                            has_acid_activation = True
                            break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_late_cyanation
        and has_amide_couplings >= 2
        and (has_nitro_reduction or has_acid_activation)
    )

    print(f"Late-stage cyanation strategy detected: {strategy_present}")
    print(f"- Late-stage cyanation: {has_late_cyanation}")
    print(f"- Amide couplings: {has_amide_couplings}")
    print(f"- Nitro reduction: {has_nitro_reduction}")
    print(f"- Acid activation: {has_acid_activation}")

    return strategy_present
