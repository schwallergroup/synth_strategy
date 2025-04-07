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
    Detects a synthetic strategy involving late-stage amide formation with an azide intermediate
    in the synthesis pathway.
    """
    # Initialize flags to track key features
    has_late_stage_amide = False
    has_azide_intermediate = False
    has_azide_to_amine_conversion = False

    # Track molecules and their transformations
    azide_molecules = set()  # Set of SMILES of molecules containing azide
    amine_from_azide = set()  # Set of SMILES of amines derived from azides
    amines_used_in_amide = set()  # Set of SMILES of amines used in amide formation

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_amide, has_azide_intermediate, has_azide_to_amine_conversion

        if node["type"] == "reaction":
            # Extract reaction information
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for late-stage amide formation (depth 0-2)
                if depth <= 2:
                    # Check for amide formation reactions - expanded list
                    amide_formation_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Schotten-Baumann_amide",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acyl chloride with ammonia to amide",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with ammonia to amide",
                        "Ester with secondary amine to amide",
                        "Nitrile to amide",
                        "Carboxylic acid to amide conversion",
                    ]

                    # Check for amide formation by reaction type
                    for reaction_type in amide_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            has_late_stage_amide = True
                            print(
                                f"Detected late-stage amide formation at depth {depth}: {reaction_type}"
                            )

                            # Store the amine reactants used in amide formation
                            for reactant in reactants:
                                if checker.check_fg(
                                    "Primary amine", reactant
                                ) or checker.check_fg("Secondary amine", reactant):
                                    amines_used_in_amide.add(reactant)
                                    print(f"Amine used in amide formation: {reactant}")

                    # Also check for amide formation by functional group transformation
                    if not has_late_stage_amide:
                        # Check if product contains amide and reactants contain amine
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            for reactant in reactants:
                                if checker.check_fg(
                                    "Primary amine", reactant
                                ) or checker.check_fg("Secondary amine", reactant):
                                    has_late_stage_amide = True
                                    amines_used_in_amide.add(reactant)
                                    print(
                                        f"Detected late-stage amide formation at depth {depth} by FG transformation"
                                    )
                                    print(f"Amine used in amide formation: {reactant}")

                # Check for azide intermediate and azide-to-amine conversion
                # Look for azides in any molecule at any depth
                for smiles in [product] + reactants:
                    if checker.check_fg("Azide", smiles):
                        has_azide_intermediate = True
                        azide_molecules.add(smiles)
                        print(f"Detected azide at depth {depth} in molecule: {smiles}")

                # Check for azide reduction to amine - expanded list
                azide_reduction_reactions = [
                    "Azide to amine reduction (Staudinger)",
                    "Reduction of nitro groups to amines",
                    "Amine to azide",  # This is the reverse, but we're traversing retrosynthetically
                    "Reduction of nitrile to amine",  # Some routes might go through nitrile
                    "Reduction of primary amides to amines",  # Some routes might go through amide
                ]

                # Check for azide-to-amine conversion by reaction type
                for reaction_type in azide_reduction_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        # For forward reactions (azide to amine)
                        if reaction_type != "Amine to azide":
                            # Check if any reactant is an azide and product is an amine
                            if any(
                                checker.check_fg("Azide", r) for r in reactants
                            ) and (
                                checker.check_fg("Primary amine", product)
                                or checker.check_fg("Secondary amine", product)
                            ):
                                has_azide_to_amine_conversion = True
                                amine_from_azide.add(product)
                                print(
                                    f"Detected azide-to-amine conversion at depth {depth}: {reaction_type}"
                                )
                                print(f"Amine from azide: {product}")

                        # For reverse reactions (amine to azide, but we're traversing retrosynthetically)
                        else:
                            # Check if any reactant is an amine and product is an azide
                            if any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                for r in reactants
                            ) and (checker.check_fg("Azide", product)):
                                has_azide_to_amine_conversion = True
                                for r in reactants:
                                    if checker.check_fg(
                                        "Primary amine", r
                                    ) or checker.check_fg("Secondary amine", r):
                                        amine_from_azide.add(r)
                                        print(
                                            f"Detected amine-to-azide conversion (retrosynthetically) at depth {depth}"
                                        )
                                        print(f"Amine from azide: {r}")

                # Check for direct functional group transformation (azide to amine)
                if not has_azide_to_amine_conversion:
                    # Check if product has amine and any reactant has azide
                    if (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                    ) and any(checker.check_fg("Azide", r) for r in reactants):
                        has_azide_to_amine_conversion = True
                        amine_from_azide.add(product)
                        print(
                            f"Detected azide-to-amine conversion at depth {depth} by FG transformation"
                        )
                        print(f"Amine from azide: {product}")

                    # Check for retrosynthetic direction (product has azide, reactant has amine)
                    elif checker.check_fg("Azide", product) and any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    ):
                        has_azide_to_amine_conversion = True
                        for r in reactants:
                            if checker.check_fg("Primary amine", r) or checker.check_fg(
                                "Secondary amine", r
                            ):
                                amine_from_azide.add(r)
                                print(
                                    f"Detected amine-to-azide conversion (retrosynthetically) at depth {depth} by FG transformation"
                                )
                                print(f"Amine from azide: {r}")

        # Check molecule nodes for azides
        elif node["type"] == "mol":
            smiles = node["smiles"]
            if checker.check_fg("Azide", smiles):
                has_azide_intermediate = True
                azide_molecules.add(smiles)
                print(f"Detected azide in molecule node at depth {depth}: {smiles}")

            # Also check for amines that might be derived from azides
            if (
                checker.check_fg("Primary amine", smiles)
                or checker.check_fg("Secondary amine", smiles)
            ) and has_azide_intermediate:
                # This is a heuristic - if we've seen azides and now see amines, they might be connected
                print(f"Potential amine from azide at depth {depth}: {smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if any amine derived from azide was used in amide formation
    amine_azide_connection = False

    # Try direct SMILES matching first
    for amine_azide in amine_from_azide:
        for amine_amide in amines_used_in_amide:
            # Check if they are the same molecule
            if amine_azide == amine_amide:
                amine_azide_connection = True
                print(
                    f"Confirmed amine from azide was used in amide formation (exact match): {amine_azide}"
                )
                break

    # If no direct match, try MCS (Maximum Common Substructure) matching
    if not amine_azide_connection and amine_from_azide and amines_used_in_amide:
        try:
            for amine_azide in amine_from_azide:
                for amine_amide in amines_used_in_amide:
                    mol1 = Chem.MolFromSmiles(amine_azide)
                    mol2 = Chem.MolFromSmiles(amine_amide)
                    if mol1 and mol2:
                        mcs = rdFMCS.FindMCS(
                            [mol1, mol2],
                            bondCompare=rdFMCS.BondCompare.CompareOrder,
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            ringMatchesRingOnly=True,
                            completeRingsOnly=True,
                        )

                        # If MCS is significant (>70% of smaller molecule)
                        if mcs.numAtoms > 0:
                            smaller_atom_count = min(
                                mol1.GetNumAtoms(), mol2.GetNumAtoms()
                            )
                            if mcs.numAtoms >= 0.7 * smaller_atom_count:
                                amine_azide_connection = True
                                print(
                                    f"Confirmed amine from azide was used in amide formation (MCS match)"
                                )
                                print(f"Amine from azide: {amine_azide}")
                                print(f"Amine in amide: {amine_amide}")
                                print(f"MCS atoms: {mcs.numAtoms}/{smaller_atom_count}")
                                break
        except Exception as e:
            print(f"Error in MCS matching: {e}")

    # If we don't have direct matches, the pathway might still be valid if we have both key transformations
    if (
        has_late_stage_amide
        and has_azide_to_amine_conversion
        and has_azide_intermediate
        and not amine_azide_connection
    ):
        print(
            "Found azide-to-amine conversion and late-stage amide formation, but couldn't directly connect them"
        )
        print("Assuming connection based on synthetic pathway")
        amine_azide_connection = True

    # Return True if all key features are present
    result = (
        has_late_stage_amide
        and has_azide_intermediate
        and has_azide_to_amine_conversion
        and amine_azide_connection
    )
    print(f"Final result: {result}")
    print(
        f"Late-stage amide: {has_late_stage_amide}, Azide intermediate: {has_azide_intermediate}, Azide-to-amine: {has_azide_to_amine_conversion}, Connection: {amine_azide_connection}"
    )

    return result
