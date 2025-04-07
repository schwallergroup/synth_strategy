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
    This function detects if the synthetic route employs a late-stage nitro reduction strategy,
    where a nitro group is carried through most of the synthesis and reduced to an amine
    only in the final step (depth 1).
    """
    # Flag to track if we found a nitro reduction at depth 1
    found_late_stage_nitro_reduction = False
    # Track if we found nitro groups in earlier steps
    found_nitro_in_earlier_steps = False
    # Track nitro-containing molecules that are later reduced
    nitro_molecules = set()

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_nitro_reduction, found_nitro_in_earlier_steps

        # Get depth from metadata if available
        if "depth" in node.get("metadata", {}):
            depth = node["metadata"]["depth"]

        print(f"Traversing node type: {node.get('type')} at depth: {depth}")

        # Check if this is a molecule node
        if node.get("type") == "mol":
            mol_smiles = node.get("smiles", "")

            # Check for nitro groups in earlier steps (depth > 1)
            if depth > 1 and checker.check_fg("Nitro group", mol_smiles):
                print(f"Found nitro group at depth {depth}: {mol_smiles}")
                found_nitro_in_earlier_steps = True
                nitro_molecules.add(mol_smiles)

        # Check if this is a reaction node
        elif node.get("type") == "reaction":
            # Check if this is the final reaction (depth 1)
            if depth == 1:
                # Extract reactants and product
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    try:
                        parts = rsmi.split(">")
                        if len(parts) >= 3:
                            reactants_smiles = parts[0].split(".")
                            product_smiles = parts[2]

                            print(f"Analyzing reaction at depth 1: {rsmi}")

                            # Check for nitro groups in reactants
                            nitro_reactants = [
                                r for r in reactants_smiles if checker.check_fg("Nitro group", r)
                            ]

                            # Check for amine groups in product
                            has_amine_product = (
                                checker.check_fg("Primary amine", product_smiles)
                                or checker.check_fg("Secondary amine", product_smiles)
                                or checker.check_fg("Tertiary amine", product_smiles)
                                or checker.check_fg("Aniline", product_smiles)
                            )

                            # Check if this is a nitro reduction reaction
                            is_nitro_reduction = checker.check_reaction(
                                "Reduction of nitro groups to amines", rsmi
                            )

                            # If not detected by reaction type, check by functional group changes
                            if not is_nitro_reduction and nitro_reactants and has_amine_product:
                                # Check if nitro group is gone in product
                                if not checker.check_fg("Nitro group", product_smiles):
                                    print(
                                        "Detected nitro reduction based on functional group changes"
                                    )
                                    is_nitro_reduction = True

                            if is_nitro_reduction and nitro_reactants:
                                print("Confirmed nitro reduction at final step")
                                # Check if any of the nitro-containing reactants were seen in earlier steps
                                for reactant in nitro_reactants:
                                    if reactant in nitro_molecules:
                                        print(
                                            f"Nitro molecule from earlier step is being reduced: {reactant}"
                                        )
                                        found_late_stage_nitro_reduction = True
                                        break

                                # If we didn't find a match in our tracked molecules, still mark as true
                                # if the reaction is definitely a nitro reduction
                                if not found_late_stage_nitro_reduction and is_nitro_reduction:
                                    print(
                                        "Nitro reduction detected but molecule wasn't tracked from earlier steps"
                                    )
                                    found_late_stage_nitro_reduction = True
                    except Exception as e:
                        print(f"Error processing reaction SMILES: {e}")

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True only if we found both a late-stage nitro reduction and nitro groups in earlier steps
    result = found_late_stage_nitro_reduction and found_nitro_in_earlier_steps
    print(
        f"Late stage nitro reduction: {found_late_stage_nitro_reduction}, Nitro in earlier steps: {found_nitro_in_earlier_steps}"
    )
    return result
