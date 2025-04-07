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
    This function detects if the synthetic route involves late-stage cyanation (nitrile introduction).
    """
    late_stage_cyanation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cyanation

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Consider reactions at depth 0, 1, or 2 (late stage)
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {depth}: {rsmi}")

                    # Check if product contains nitrile
                    product_has_nitrile = checker.check_fg("Nitrile", product_smiles)
                    print(f"Product has nitrile: {product_has_nitrile}")

                    if product_has_nitrile:
                        # Check if any reactant has a true nitrile functional group (not a cyanide reagent)
                        reactants_with_nitrile = []
                        for r_smi in reactants_smiles:
                            if checker.check_fg("Nitrile", r_smi):
                                # Check if this is a cyanide reagent rather than a molecule with nitrile group
                                mol = Chem.MolFromSmiles(r_smi)
                                if mol:
                                    # Common cyanide reagents like [Cu]C#N, NaC#N, KC#N, etc.
                                    if "[Cu]" in r_smi and "#" in r_smi:
                                        print(
                                            f"Detected copper cyanide reagent: {r_smi}"
                                        )
                                        continue
                                    if (
                                        "Na" in r_smi
                                        and "#" in r_smi
                                        and len(r_smi) < 10
                                    ):
                                        print(
                                            f"Detected sodium cyanide reagent: {r_smi}"
                                        )
                                        continue
                                    if (
                                        "K" in r_smi
                                        and "#" in r_smi
                                        and len(r_smi) < 10
                                    ):
                                        print(
                                            f"Detected potassium cyanide reagent: {r_smi}"
                                        )
                                        continue
                                    if "TMS" in r_smi and "#" in r_smi:
                                        print(f"Detected TMSCN reagent: {r_smi}")
                                        continue

                                    # If we get here, it's likely a true nitrile-containing molecule
                                    reactants_with_nitrile.append(r_smi)
                                    print(
                                        f"Reactant has nitrile functional group: {r_smi}"
                                    )

                        # Check for cyanide reagents explicitly
                        has_cyanide_reagent = False
                        for r_smi in reactants_smiles:
                            # Common cyanide reagents
                            if (
                                "[Cu]C#N" in r_smi
                                or "[Cu][C]#[N]" in r_smi
                                or "NaCN" in r_smi
                                or "KCN" in r_smi
                                or "TMSCN" in r_smi
                                or ("[Cu]" in r_smi and "#" in r_smi)
                            ):
                                has_cyanide_reagent = True
                                print(f"Detected cyanide reagent: {r_smi}")
                                break

                        # If product has nitrile, no reactants have true nitrile groups, and there's a cyanide reagent
                        # OR if it's a known cyanation reaction, then it's a cyanation reaction
                        if (
                            product_has_nitrile
                            and len(reactants_with_nitrile) == 0
                            and has_cyanide_reagent
                        ):
                            print(
                                f"Late-stage cyanation detected in reaction at depth {depth}: {rsmi}"
                            )
                            print(f"Product: {product_smiles}")
                            print(f"Reactants: {reactants_smiles}")
                            late_stage_cyanation = True
                            return

                        # Check if this is a known reaction involving nitriles
                        nitrile_related_reactions = [
                            "Aryl halide to carboxylic acid",  # Can involve nitrile intermediate
                            "Oxidation of nitrile to carboxylic acid",
                            "Nitrile to amide",
                            "Reduction of nitrile to amine",
                        ]

                        for rxn_type in nitrile_related_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Detected nitrile-related reaction: {rxn_type}")
                                # These are transformations OF nitriles, not cyanation reactions
                                break

                        # Check if reactants have functional groups commonly converted to nitriles
                        # and if there's no true nitrile in reactants but there is in product
                        if len(reactants_with_nitrile) == 0 and product_has_nitrile:
                            reactants_have_convertible_group = False
                            for r_smi in reactants_smiles:
                                if (
                                    checker.check_fg("Primary halide", r_smi)
                                    or checker.check_fg("Secondary halide", r_smi)
                                    or checker.check_fg("Tertiary halide", r_smi)
                                    or checker.check_fg("Aromatic halide", r_smi)
                                    or checker.check_fg("Aldehyde", r_smi)
                                ):
                                    reactants_have_convertible_group = True
                                    print(f"Reactant has convertible group: {r_smi}")
                                    break

                            if reactants_have_convertible_group:
                                print(
                                    f"Late-stage cyanation detected in reaction at depth {depth}: {rsmi}"
                                )
                                print(f"Product: {product_smiles}")
                                print(f"Reactants: {reactants_smiles}")
                                late_stage_cyanation = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    print("Starting late-stage cyanation detection...")
    dfs_traverse(route)
    print(f"Late-stage cyanation detected: {late_stage_cyanation}")
    return late_stage_cyanation
