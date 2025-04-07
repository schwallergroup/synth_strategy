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
    This function detects a homologation strategy involving nitrile installation,
    followed by hydrolysis to carboxylic acid, and subsequent transformations.
    """
    # Track if we found each step in the sequence
    nitrile_installation = {"found": False, "depth": -1, "atom_indices": None}
    nitrile_hydrolysis = {"found": False, "depth": -1, "atom_indices": None}
    acid_transformation = {"found": False, "depth": -1, "atom_indices": None}

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_installation, nitrile_hydrolysis, acid_transformation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")
                product = product_part

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for nitrile installation (mesylate/other leaving group → nitrile)
                if not nitrile_installation["found"]:
                    # Check if product has nitrile
                    if checker.check_fg("Nitrile", product):
                        print(f"Product has nitrile: {product}")

                        # Get nitrile atom indices in product
                        nitrile_indices = checker.get_fg_atom_indices("Nitrile", product)

                        # Check if reactants have leaving groups
                        has_leaving_group = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Mesylate", reactant)
                                or checker.check_fg("Primary halide", reactant)
                                or checker.check_fg("Secondary halide", reactant)
                                or checker.check_fg("Tertiary halide", reactant)
                                or checker.check_fg("Triflate", reactant)
                            ):
                                has_leaving_group = True
                                print(f"Found leaving group in reactant: {reactant}")
                                break

                        if has_leaving_group:
                            nitrile_installation["found"] = True
                            nitrile_installation["depth"] = depth
                            nitrile_installation["atom_indices"] = nitrile_indices
                            print(f"Found nitrile installation step at depth {depth}")

                # Check for nitrile hydrolysis (nitrile → carboxylic acid)
                if not nitrile_hydrolysis["found"]:
                    # Check if reactant has nitrile
                    nitrile_reactant = None
                    for reactant in reactants:
                        if checker.check_fg("Nitrile", reactant):
                            nitrile_reactant = reactant
                            print(f"Reactant has nitrile: {reactant}")
                            break

                    # Check if product has carboxylic acid
                    if nitrile_reactant and checker.check_fg("Carboxylic acid", product):
                        print(f"Product has carboxylic acid: {product}")

                        # Get carboxylic acid atom indices in product
                        acid_indices = checker.get_fg_atom_indices("Carboxylic acid", product)

                        nitrile_hydrolysis["found"] = True
                        nitrile_hydrolysis["depth"] = depth
                        nitrile_hydrolysis["atom_indices"] = acid_indices
                        print(f"Found nitrile hydrolysis step at depth {depth}")

                # Check for subsequent transformation of carboxylic acid
                if not acid_transformation["found"]:
                    # Check if reactant has carboxylic acid
                    acid_reactant = None
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            acid_reactant = reactant
                            print(f"Reactant has carboxylic acid: {reactant}")
                            break

                    if acid_reactant:
                        # Check for common acid transformations
                        transformation_found = False

                        # Check for specific reaction types
                        if (
                            checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                            or checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                            )
                            or checker.check_reaction(
                                "Reduction of carboxylic acid to primary alcohol", rsmi
                            )
                            or checker.check_reaction("Decarboxylation", rsmi)
                        ):
                            transformation_found = True

                        # Also check for product functional groups that indicate transformation
                        elif (
                            checker.check_fg("Ester", product)
                            or checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                            or checker.check_fg("Primary alcohol", product)
                        ):
                            transformation_found = True

                        if transformation_found:
                            acid_transformation["found"] = True
                            acid_transformation["depth"] = depth
                            print(f"Found carboxylic acid transformation step at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if all steps were found
    all_steps_found = (
        nitrile_installation["found"]
        and nitrile_hydrolysis["found"]
        and acid_transformation["found"]
    )

    # Check if steps are in the correct sequence (in retrosynthetic order)
    correct_sequence = False
    if all_steps_found:
        # In retrosynthetic order, acid transformation should be at lowest depth,
        # followed by nitrile hydrolysis, then nitrile installation at highest depth
        correct_sequence = (
            acid_transformation["depth"]
            < nitrile_hydrolysis["depth"]
            < nitrile_installation["depth"]
        )

    strategy_detected = all_steps_found and correct_sequence

    print(f"Nitrile homologation strategy detected: {strategy_detected}")
    print(
        f"Steps found: nitrile installation: {nitrile_installation['found']} (depth {nitrile_installation['depth']}), "
        f"nitrile hydrolysis: {nitrile_hydrolysis['found']} (depth {nitrile_hydrolysis['depth']}), "
        f"acid transformation: {acid_transformation['found']} (depth {acid_transformation['depth']})"
    )
    print(f"Correct sequence: {correct_sequence}")

    return strategy_detected
