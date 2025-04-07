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
    Detects a synthetic strategy involving late-stage amide coupling with preceding
    functional group activation (acid → acid chloride → amide) and protection/deprotection.
    """
    # Track if we found the key elements of the strategy
    found_amide_coupling = False
    found_acid_chloride_formation = False
    found_ester_hydrolysis = False
    found_acid_activation = False  # Track any acid activation method

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_coupling, found_acid_chloride_formation, found_ester_hydrolysis, found_acid_activation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = reactants_str.split(".")

                # Check for amide coupling at late stage (depth 0-2)
                if depth <= 2:
                    # Check if this is a Schotten-Baumann amide formation or similar reaction
                    if (
                        checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                            rsmi,
                        )
                        or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            rsmi,
                        )
                    ):
                        found_amide_coupling = True
                        print(f"Found late-stage amide coupling at depth {depth}")
                    else:
                        # Fallback to functional group checking
                        has_acid_halide = False
                        has_carboxylic_acid = False
                        has_amine = False
                        for reactant in reactants:
                            if checker.check_fg("Acyl halide", reactant):
                                has_acid_halide = True
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_carboxylic_acid = True
                            if checker.check_fg(
                                "Primary amine", reactant
                            ) or checker.check_fg("Secondary amine", reactant):
                                has_amine = True

                        has_amide_product = (
                            checker.check_fg("Primary amide", product_str)
                            or checker.check_fg("Secondary amide", product_str)
                            or checker.check_fg("Tertiary amide", product_str)
                        )

                        if (
                            (has_acid_halide or has_carboxylic_acid)
                            and has_amine
                            and has_amide_product
                        ):
                            found_amide_coupling = True
                            print(
                                f"Found late-stage amide coupling (FG check) at depth {depth}"
                            )

                # Check for acid chloride formation or other activation (depth 1-3)
                if 1 <= depth <= 3:
                    # Check for specific reaction types
                    if checker.check_reaction(
                        "Acyl chlorides from alcohols", rsmi
                    ) or checker.check_reaction("Alcohol to chloride_SOCl2", rsmi):
                        found_acid_chloride_formation = True
                        found_acid_activation = True
                        print(f"Found acid chloride formation at depth {depth}")
                    else:
                        # Fallback to functional group checking
                        has_carboxylic_acid = False
                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_carboxylic_acid = True

                        has_acid_chloride_product = checker.check_fg(
                            "Acyl halide", product_str
                        )

                        if has_carboxylic_acid and has_acid_chloride_product:
                            found_acid_chloride_formation = True
                            found_acid_activation = True
                            print(
                                f"Found acid chloride formation (FG check) at depth {depth}"
                            )

                # Check for ester hydrolysis (deprotection) or other acid generation (depth 2-4)
                if 2 <= depth <= 4:
                    # Check for specific reaction types
                    if (
                        checker.check_reaction(
                            "Ester saponification (alkyl deprotection)", rsmi
                        )
                        or checker.check_reaction(
                            "Ester saponification (methyl deprotection)", rsmi
                        )
                        or checker.check_reaction(
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                            rsmi,
                        )
                        or checker.check_reaction("COOH ethyl deprotection", rsmi)
                        or checker.check_reaction(
                            "Deprotection of carboxylic acid", rsmi
                        )
                        or checker.check_reaction(
                            "Oxidation of aldehydes to carboxylic acids", rsmi
                        )
                        or checker.check_reaction(
                            "Oxidation of nitrile to carboxylic acid", rsmi
                        )
                    ):
                        found_ester_hydrolysis = True
                        print(f"Found carboxylic acid generation at depth {depth}")
                    else:
                        # Fallback to functional group checking
                        has_protected_acid = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Ester", reactant)
                                or checker.check_fg("Nitrile", reactant)
                                or checker.check_fg("Aldehyde", reactant)
                            ):
                                has_protected_acid = True

                        has_carboxylic_acid_product = checker.check_fg(
                            "Carboxylic acid", product_str
                        )

                        if has_protected_acid and has_carboxylic_acid_product:
                            found_ester_hydrolysis = True
                            print(
                                f"Found carboxylic acid generation (FG check) at depth {depth}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Print overall results for debugging
    print(f"Amide coupling found: {found_amide_coupling}")
    print(f"Acid chloride formation found: {found_acid_chloride_formation}")
    print(f"Ester hydrolysis/acid generation found: {found_ester_hydrolysis}")

    # Return True if we found amide coupling and either acid activation or ester hydrolysis
    # This captures the essence of the strategy while allowing for variations
    return found_amide_coupling and (
        found_acid_chloride_formation or found_ester_hydrolysis
    )
