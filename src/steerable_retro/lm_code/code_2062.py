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
    This function detects a sequence of nitrogen oxidation state changes:
    nitro → amine → hydrazine → heterocycle (pyrazole)

    In retrosynthetic traversal, we expect to find these in reverse order:
    pyrazole (target) → hydrazine → amine → nitro (starting material)
    """
    # Track the stages we've found
    nitro_found = False
    amine_found = False
    hydrazine_found = False
    pyrazole_found = False

    # Track the order (depth) of transformations
    # In retrosynthesis, lower depth = later stage (closer to target)
    nitro_depth = -1
    amine_depth = -1
    hydrazine_depth = -1
    pyrazole_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal nitro_found, amine_found, hydrazine_found, pyrazole_found
        nonlocal nitro_depth, amine_depth, hydrazine_depth, pyrazole_depth

        if node["type"] == "mol":
            # Check for functional groups in molecules
            mol_smiles = node["smiles"]

            # Record the occurrence of each group, updating if found at an earlier stage
            if checker.check_fg("Nitro group", mol_smiles) and (
                nitro_depth == -1 or depth < nitro_depth
            ):
                print(f"Found nitro group in molecule at depth {depth}")
                nitro_found = True
                nitro_depth = depth

            if checker.check_fg("Primary amine", mol_smiles) and (
                amine_depth == -1 or depth < amine_depth
            ):
                print(f"Found primary amine in molecule at depth {depth}")
                amine_found = True
                amine_depth = depth

            if checker.check_fg("Hydrazine", mol_smiles) and (
                hydrazine_depth == -1 or depth < hydrazine_depth
            ):
                print(f"Found hydrazine in molecule at depth {depth}")
                hydrazine_found = True
                hydrazine_depth = depth

            if checker.check_ring("pyrazole", mol_smiles) and (
                pyrazole_depth == -1 or depth < pyrazole_depth
            ):
                print(f"Found pyrazole in molecule at depth {depth}")
                pyrazole_found = True
                pyrazole_depth = depth

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction to amine
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction to amine at depth {depth}")
                    if any(
                        checker.check_fg("Nitro group", r) for r in reactants
                    ) and checker.check_fg("Primary amine", product):
                        if amine_depth == -1 or depth < amine_depth:
                            amine_found = True
                            amine_depth = depth

                # Check for amine to hydrazine transformation
                if checker.check_reaction("Hydrazine synthesis from amine", rsmi):
                    print(
                        f"Found amine to hydrazine transformation (specific reaction) at depth {depth}"
                    )
                    if hydrazine_depth == -1 or depth < hydrazine_depth:
                        hydrazine_found = True
                        hydrazine_depth = depth

                # Generic check for amine to hydrazine transformation
                if (
                    any(checker.check_fg("Primary amine", r) for r in reactants)
                    and checker.check_fg("Hydrazine", product)
                    and not any(checker.check_fg("Hydrazine", r) for r in reactants)
                ):
                    print(f"Found amine to hydrazine transformation at depth {depth}")
                    if hydrazine_depth == -1 or depth < hydrazine_depth:
                        hydrazine_found = True
                        hydrazine_depth = depth

                # Check for hydrazine to pyrazole transformation
                if checker.check_reaction("pyrazole", rsmi) or checker.check_reaction(
                    "Pyrazole formation", rsmi
                ):
                    print(f"Found pyrazole formation reaction at depth {depth}")
                    if any(
                        checker.check_fg("Hydrazine", r) for r in reactants
                    ) and checker.check_ring("pyrazole", product):
                        if pyrazole_depth == -1 or depth < pyrazole_depth:
                            pyrazole_found = True
                            pyrazole_depth = depth

                # Generic check for hydrazine to pyrazole transformation
                if (
                    any(checker.check_fg("Hydrazine", r) for r in reactants)
                    and checker.check_ring("pyrazole", product)
                    and not any(checker.check_ring("pyrazole", r) for r in reactants)
                ):
                    print(f"Found hydrazine to pyrazole transformation at depth {depth}")
                    if pyrazole_depth == -1 or depth < pyrazole_depth:
                        pyrazole_found = True
                        pyrazole_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro found: {nitro_found} at depth {nitro_depth}")
    print(f"Amine found: {amine_found} at depth {amine_depth}")
    print(f"Hydrazine found: {hydrazine_found} at depth {hydrazine_depth}")
    print(f"Pyrazole found: {pyrazole_found} at depth {pyrazole_depth}")

    # Check if we found all stages in the correct order for retrosynthesis
    # In retrosynthetic traversal, target (pyrazole) has lowest depth, starting material (nitro) has highest
    # Allow for equal depths in case multiple functional groups are in the same molecule
    correct_sequence = (
        nitro_found
        and amine_found
        and hydrazine_found
        and pyrazole_found
        and pyrazole_depth <= hydrazine_depth <= amine_depth <= nitro_depth
    )

    print(f"Correct sequence: {correct_sequence}")
    return correct_sequence
