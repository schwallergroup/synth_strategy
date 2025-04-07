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
    Detects if the synthesis uses late-stage functionalization, where significant
    functional groups or structural elements are introduced in the final steps.
    """
    late_stage_depth_threshold = 1  # Consider depth 0-1 as "late stage"

    late_stage_functionalization_detected = False

    def dfs_traverse(node):
        nonlocal late_stage_functionalization_detected

        if node["type"] == "reaction":
            depth = node.get("metadata", {}).get("depth", -1)

            if depth <= late_stage_depth_threshold:
                try:
                    rsmi = node["metadata"].get("rsmi", "")
                    if not rsmi:
                        return

                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    if not reactants_part or not product_part:
                        return

                    reactants = reactants_part.split(".")

                    # Check for important functional groups
                    significant_fgs = [
                        "Tetrazole",
                        "Triazole",
                        "Sulfonamide",
                        "Nitro group",
                        "Azide",
                        "Isocyanate",
                        "Isothiocyanate",
                        "Nitrile",
                        "Carboxylic acid",
                        "Ester",
                        "Amide",
                        "Urea",
                        "Thiourea",
                    ]

                    for fg_name in significant_fgs:
                        if checker.check_fg(fg_name, product_part):
                            fg_in_all_reactants = all(
                                checker.check_fg(fg_name, r) for r in reactants
                            )
                            if not fg_in_all_reactants:
                                late_stage_functionalization_detected = True
                                print(
                                    f"Detected late-stage {fg_name} introduction at depth {depth}"
                                )

                    # Check for important ring structures
                    important_rings = [
                        "tetrazole",
                        "triazole",
                        "isoxazole",
                        "oxadiazole",
                        "thiazole",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "pyrimidine",
                        "pyridine",
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                    ]

                    for ring_name in important_rings:
                        if checker.check_ring(ring_name, product_part):
                            ring_in_all_reactants = all(
                                checker.check_ring(ring_name, r) for r in reactants
                            )
                            if not ring_in_all_reactants:
                                late_stage_functionalization_detected = True
                                print(
                                    f"Detected late-stage {ring_name} ring introduction at depth {depth}"
                                )

                    # Check for specific late-stage functionalization reactions
                    late_stage_rxn_types = [
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        "tetrazole_terminal",
                        "tetrazole_connect_regioisomere_1",
                        "tetrazole_connect_regioisomere_2",
                        "Azide-nitrile click cycloaddition to tetrazole",
                        "Azide-nitrile click cycloaddition to triazole",
                        "pyrazole",
                        "oxadiazole",
                        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                        "Urea synthesis via isocyanate and primary amine",
                        "Urea synthesis via isocyanate and secondary amine",
                    ]

                    for rxn_type in late_stage_rxn_types:
                        if checker.check_reaction(rxn_type, rsmi):
                            late_stage_functionalization_detected = True
                            print(f"Detected late-stage reaction: {rxn_type} at depth {depth}")

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return late_stage_functionalization_detected
