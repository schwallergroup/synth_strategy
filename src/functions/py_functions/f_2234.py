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
    Detects if the synthesis uses a protection-deprotection sequence
    with Cbz or Boc protecting groups
    """
    protection_depths = []
    deprotection_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal protection_depths, deprotection_depths

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection reactions (Boc)
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction(
                        "Boc amine protection with Boc anhydride", rsmi
                    )
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction(
                        "Boc amine protection of secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Boc amine protection of primary amine", rsmi
                    )
                ):
                    print(f"Found Boc protection step at depth {depth}, rsmi: {rsmi}")
                    protection_depths.append(depth)

                # Check for deprotection reactions (Boc)
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction(
                        "Boc amine deprotection of guanidine", rsmi
                    )
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    print(f"Found Boc deprotection step at depth {depth}, rsmi: {rsmi}")
                    deprotection_depths.append(depth)

                # Check for Cbz protection
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol is not None:
                    has_cbz_in_product = checker.check_fg("Carbamic ester", product)

                    if has_cbz_in_product:
                        # Check if reactants don't have Cbz
                        all_reactants_without_cbz = True
                        for reactant in reactants:
                            if checker.check_fg("Carbamic ester", reactant):
                                all_reactants_without_cbz = False
                                break

                        if all_reactants_without_cbz:
                            # Check if the reaction involves an amine
                            for reactant in reactants:
                                if (
                                    checker.check_fg("Primary amine", reactant)
                                    or checker.check_fg("Secondary amine", reactant)
                                    or checker.check_fg("Aniline", reactant)
                                ):
                                    # Check for benzyl chloroformate or similar reagents
                                    for r in reactants:
                                        if (
                                            "ClC(=O)Oc1ccccc1" in r
                                            or "ClCOOc1ccccc1" in r
                                            or "C(=O)Cl" in r
                                        ):
                                            print(
                                                f"Found Cbz protection step at depth {depth}, rsmi: {rsmi}"
                                            )
                                            protection_depths.append(depth)
                                            break

                                    # If we didn't find a specific reagent but the pattern matches
                                    if depth not in protection_depths:
                                        print(
                                            f"Found Cbz protection step at depth {depth}, rsmi: {rsmi}"
                                        )
                                        protection_depths.append(depth)
                                    break

                # Check for Cbz deprotection
                for reactant in reactants:
                    if checker.check_fg("Carbamic ester", reactant):
                        # Check if product doesn't have Cbz
                        if prod_mol is not None and not checker.check_fg(
                            "Carbamic ester", product
                        ):
                            # Check for specific deprotection reactions
                            if (
                                checker.check_reaction(
                                    "Hydrogenolysis of amides/imides/carbamates", rsmi
                                )
                                or checker.check_reaction(
                                    "Hydrolysis of amides/imides/carbamates", rsmi
                                )
                                or checker.check_reaction(
                                    "Carboxyl benzyl deprotection", rsmi
                                )
                                or checker.check_reaction(
                                    "Hydroxyl benzyl deprotection", rsmi
                                )
                            ):
                                print(
                                    f"Found Cbz deprotection step at depth {depth}, rsmi: {rsmi}"
                                )
                                deprotection_depths.append(depth)
                                break

                            # Check for general hydrogenolysis patterns
                            for r in reactants:
                                if "H2" in r or "[H][H]" in r or "Pd" in r or "Pt" in r:
                                    print(
                                        f"Found Cbz deprotection step (hydrogenolysis) at depth {depth}, rsmi: {rsmi}"
                                    )
                                    deprotection_depths.append(depth)
                                    break

                            # Check for hydrolysis patterns
                            if (
                                "OH" in rsmi
                                or "H2O" in rsmi
                                or "NaOH" in rsmi
                                or "KOH" in rsmi
                            ):
                                if depth not in deprotection_depths:
                                    print(
                                        f"Found Cbz deprotection step (hydrolysis) at depth {depth}, rsmi: {rsmi}"
                                    )
                                    deprotection_depths.append(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if we have both protection and deprotection
    has_protection = len(protection_depths) > 0
    has_deprotection = len(deprotection_depths) > 0

    # Check if there's a valid sequence (protection should come before deprotection in the synthesis)
    valid_sequence = False
    if has_protection and has_deprotection:
        # In retrosynthesis, higher depth means earlier in the actual synthesis
        # So we need to find if any protection depth is higher than any deprotection depth
        for p_depth in protection_depths:
            for d_depth in deprotection_depths:
                if (
                    p_depth > d_depth
                ):  # Protection occurs earlier in synthesis than deprotection
                    valid_sequence = True
                    break
            if valid_sequence:
                break

    print(
        f"Protection found: {has_protection}, Deprotection found: {has_deprotection}, Valid sequence: {valid_sequence}"
    )
    print(
        f"Protection depths: {protection_depths}, Deprotection depths: {deprotection_depths}"
    )

    # Return True if we have both protection and deprotection in the correct order
    return has_protection and has_deprotection and valid_sequence
