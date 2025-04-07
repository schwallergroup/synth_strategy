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
    This function detects if the synthetic route involves multiple SNAr reactions (aromatic halide with amine).
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            print(f"Analyzing reaction: {rsmi}")

            # Check if this is a nucleophilic aromatic substitution reaction
            is_snar = False

            # First check if the reaction matches known SNAr patterns
            snar_reaction_types = [
                "Ullmann-Goldberg Substitution amine",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "heteroaromatic_nuc_sub",
                "nucl_sub_aromatic_ortho_nitro",
                "nucl_sub_aromatic_para_nitro",
                "Goldberg coupling",
                "Goldberg coupling aryl amine-aryl chloride",
                "Goldberg coupling aryl amide-aryl chloride",
                "Ullmann condensation",
                "N-arylation_heterocycles",
                "Buchwald-Hartwig",
            ]

            for reaction_type in snar_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    is_snar = True
                    print(f"Detected SNAr reaction via reaction pattern: {reaction_type}")
                    break

            # If not identified by reaction pattern, check for functional group transformations
            if not is_snar:
                # Check for direct evidence of SNAr by looking for aromatic C-N bond formation
                has_aromatic_halide = False
                has_amine = False
                has_activating_group = False

                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aromatic_halide = True
                        print(f"Found aromatic halide in reactant: {reactant}")

                        # Check for activating groups that facilitate SNAr
                        if (
                            checker.check_fg("Nitro group", reactant)
                            or checker.check_fg("Nitrile", reactant)
                            or checker.check_fg("Trifluoro group", reactant)
                        ):
                            has_activating_group = True
                            print(f"Found activating group in reactant: {reactant}")

                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Tertiary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                    ):
                        has_amine = True
                        print(f"Found amine in reactant: {reactant}")

                # Check if product has N-aryl bond that wasn't in reactants
                has_n_aryl_product = False
                if (
                    checker.check_fg("Aniline", product)
                    or checker.check_fg("Tertiary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_ring("pyrrole", product)
                    or checker.check_ring("indole", product)
                    or checker.check_ring("pyrazole", product)
                    or checker.check_ring("imidazole", product)
                    or checker.check_ring("triazole", product)
                ):

                    # Check if this N-aryl bond is new (not present in reactants)
                    n_aryl_in_reactants = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Aniline", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            and checker.check_ring("benzene", reactant)
                        ):
                            n_aryl_in_reactants = True

                    if not n_aryl_in_reactants:
                        has_n_aryl_product = True
                        print(f"Found new N-aryl bond in product: {product}")

                # Determine if this is an SNAr reaction based on functional group analysis
                if has_aromatic_halide and has_amine and has_n_aryl_product:
                    is_snar = True
                    print(f"Detected SNAr reaction via functional group analysis")
                    # If we have an activating group, it's even more likely to be SNAr
                    if has_activating_group:
                        print(f"Confirmed SNAr by presence of activating group")

            if is_snar:
                snar_count += 1
                print(f"SNAr reaction count: {snar_count}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total SNAr reactions found: {snar_count}")

    return snar_count >= 2
