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
    Detects a linear synthesis with sequential functionalization of a core structure.

    A sequential functionalization strategy involves a series of reactions that progressively
    modify a core structure by adding or transforming functional groups in a logical sequence.
    """
    # Track functionalization steps with their depths
    functionalization_steps = []

    def dfs_traverse(node, depth=0):
        # Process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for common functionalization reactions
            reaction_detected = False

            # Amide formation
            if (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                )
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                )
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                )
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("{Schotten-Baumann_amide}", rsmi)
            ):
                functionalization_steps.append(("amide_formation", depth))
                reaction_detected = True
                print(f"Detected amide formation at depth {depth}: {rsmi}")

            # Ester formation
            elif (
                checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                or checker.check_reaction("Transesterification", rsmi)
                or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                or checker.check_reaction("Acetic anhydride and alcohol to ester", rsmi)
            ):
                functionalization_steps.append(("ester_formation", depth))
                reaction_detected = True
                print(f"Detected ester formation at depth {depth}: {rsmi}")

            # Phenol alkylation / Ether formation
            elif (
                checker.check_reaction("Williamson Ether Synthesis", rsmi)
                or checker.check_reaction("{Williamson ether}", rsmi)
                or checker.check_reaction("Alcohol to ether", rsmi)
            ):
                # Check if one of the reactants contains a phenol
                for reactant in reactants:
                    if checker.check_fg("Phenol", reactant):
                        functionalization_steps.append(("phenol_alkylation", depth))
                        reaction_detected = True
                        print(f"Detected phenol alkylation at depth {depth}: {rsmi}")
                        break

                # If not a phenol alkylation but still an ether formation
                if not reaction_detected:
                    functionalization_steps.append(("ether_formation", depth))
                    reaction_detected = True
                    print(f"Detected ether formation at depth {depth}: {rsmi}")

            # Nucleophilic aromatic substitution
            elif (
                checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
            ):
                functionalization_steps.append(("nucleophilic_substitution", depth))
                reaction_detected = True
                print(f"Detected nucleophilic substitution at depth {depth}: {rsmi}")

            # Halogenation
            elif (
                checker.check_reaction("Aromatic chlorination", rsmi)
                or checker.check_reaction("Aromatic bromination", rsmi)
                or checker.check_reaction("Aromatic fluorination", rsmi)
                or checker.check_reaction("Aromatic iodination", rsmi)
                or checker.check_reaction("Chlorination", rsmi)
                or checker.check_reaction("Bromination", rsmi)
                or checker.check_reaction("Fluorination", rsmi)
                or checker.check_reaction("Iodination", rsmi)
            ):
                functionalization_steps.append(("halogenation", depth))
                reaction_detected = True
                print(f"Detected halogenation at depth {depth}: {rsmi}")

            # Oxidation reactions
            elif (
                checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                or checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi)
                or checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                )
                or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi)
                or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
            ):
                functionalization_steps.append(("oxidation", depth))
                reaction_detected = True
                print(f"Detected oxidation reaction at depth {depth}: {rsmi}")

            # Reduction reactions
            elif (
                checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi)
                or checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                or checker.check_reaction("Reduction of carboxylic acid to primary alcohol", rsmi)
                or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                or checker.check_reaction("Reduction of primary amides to amines", rsmi)
                or checker.check_reaction("Reduction of nitro groups to amines", rsmi)
            ):
                functionalization_steps.append(("reduction", depth))
                reaction_detected = True
                print(f"Detected reduction reaction at depth {depth}: {rsmi}")

            # Alkylation reactions
            elif (
                checker.check_reaction("Alkylation of amines", rsmi)
                or checker.check_reaction("N-alkylation of primary amines with alkyl halides", rsmi)
                or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )
                or checker.check_reaction("S-alkylation of thiols", rsmi)
                or checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rsmi
                )
            ):
                functionalization_steps.append(("alkylation", depth))
                reaction_detected = True
                print(f"Detected alkylation reaction at depth {depth}: {rsmi}")

            # Protection reactions
            elif (
                checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                or checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Protection of carboxylic acid", rsmi)
            ):
                functionalization_steps.append(("protection", depth))
                reaction_detected = True
                print(f"Detected protection reaction at depth {depth}: {rsmi}")

            # Deprotection reactions
            elif (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Deprotection of carboxylic acid", rsmi)
                or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
            ):
                functionalization_steps.append(("deprotection", depth))
                reaction_detected = True
                print(f"Detected deprotection reaction at depth {depth}: {rsmi}")

            # If no specific reaction was detected but there's a transformation of functional groups
            if not reaction_detected:
                # Check for functional group changes by comparing reactants and product
                for reactant in reactants:
                    # Check for functional groups that disappear
                    for fg in [
                        "Phenol",
                        "Carboxylic acid",
                        "Ester",
                        "Amide",
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Ether",
                        "Nitrile",
                        "Nitro group",
                        "Aldehyde",
                        "Ketone",
                        "Acyl halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                        "Aromatic halide",
                    ]:
                        if checker.check_fg(fg, reactant) and not checker.check_fg(fg, product):
                            # Check if a new functional group appears in the product
                            new_fg_found = False
                            for new_fg in [
                                "Phenol",
                                "Carboxylic acid",
                                "Ester",
                                "Amide",
                                "Primary amine",
                                "Secondary amine",
                                "Tertiary amine",
                                "Primary alcohol",
                                "Secondary alcohol",
                                "Tertiary alcohol",
                                "Ether",
                                "Nitrile",
                                "Nitro group",
                                "Aldehyde",
                                "Ketone",
                            ]:
                                if not checker.check_fg(new_fg, reactant) and checker.check_fg(
                                    new_fg, product
                                ):
                                    functionalization_steps.append(
                                        (f"{fg.lower()}_to_{new_fg.lower()}", depth)
                                    )
                                    print(
                                        f"Detected {fg.lower()} to {new_fg.lower()} transformation at depth {depth}: {rsmi}"
                                    )
                                    reaction_detected = True
                                    new_fg_found = True
                                    break

                            # If no new functional group was found but the old one disappeared
                            if not new_fg_found:
                                functionalization_steps.append(
                                    (f"{fg.lower()}_transformation", depth)
                                )
                                print(
                                    f"Detected {fg.lower()} transformation at depth {depth}: {rsmi}"
                                )
                                reaction_detected = True
                            break

                    if reaction_detected:
                        break

                # If still no reaction detected, check for new functional groups appearing
                if not reaction_detected:
                    for fg in [
                        "Phenol",
                        "Carboxylic acid",
                        "Ester",
                        "Amide",
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Ether",
                        "Nitrile",
                        "Nitro group",
                        "Aldehyde",
                        "Ketone",
                    ]:
                        if checker.check_fg(fg, product) and not any(
                            checker.check_fg(fg, r) for r in reactants
                        ):
                            functionalization_steps.append((f"{fg.lower()}_formation", depth))
                            print(f"Detected {fg.lower()} formation at depth {depth}: {rsmi}")
                            reaction_detected = True
                            break

        # Traverse children (in retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort functionalization steps by depth to check if they're sequential
    functionalization_steps.sort(key=lambda x: x[1])

    # Check if we have at least 4 different functionalization steps
    unique_steps = set(step[0] for step in functionalization_steps)

    if len(unique_steps) >= 4:
        # Check if the steps are generally sequential
        # We don't require strictly increasing depths, just a logical progression
        depths = [step[1] for step in functionalization_steps]

        # Check if depths are non-decreasing or have at most one decrease
        decreases = sum(1 for i in range(len(depths) - 1) if depths[i] > depths[i + 1])
        is_sequential = decreases <= 1

        if is_sequential:
            print(f"Found sequential functionalization with steps: {unique_steps}")
            print(f"Steps in order: {functionalization_steps}")
            return True
        else:
            print(f"Found {len(unique_steps)} functionalization steps, but they are not sequential")
    else:
        print(f"Found only {len(unique_steps)} different functionalization steps, need at least 4")

    return False
