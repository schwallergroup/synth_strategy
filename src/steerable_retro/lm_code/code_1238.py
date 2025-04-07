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
    Detects if the synthesis route uses a late-stage amide coupling strategy.
    Late-stage means the final step (depth 0 or 1) involves an amide coupling.
    """
    amide_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_found

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if this is at depth 0 or 1 (late stage)
            # Depth 1 is included because the root node might be a molecule
            if depth <= 1:
                print(f"This is a late-stage reaction (depth {depth})")

                # Check for various amide coupling reactions
                is_amide_coupling = (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                )

                print(f"Is amide coupling reaction: {is_amide_coupling}")

                # Verify product contains an amide group
                has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                print(f"Product contains amide group: {has_amide}")

                # Check reactants for carboxylic acid/derivatives and amines
                has_carboxylic_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                has_ester = any(checker.check_fg("Ester", r) for r in reactants)

                has_amine = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    or checker.check_fg("Aniline", r)
                    for r in reactants
                )

                print(f"Reactants contain carboxylic acid: {has_carboxylic_acid}")
                print(f"Reactants contain acyl halide: {has_acyl_halide}")
                print(f"Reactants contain ester: {has_ester}")
                print(f"Reactants contain amine: {has_amine}")

                # Additional check for amide formation
                if is_amide_coupling and has_amide:
                    print(f"Detected late-stage amide coupling at depth {depth}: {rsmi}")
                    amide_coupling_found = True
                    return  # Early return once found
                elif (
                    has_amide
                    and (has_carboxylic_acid or has_acyl_halide or has_ester)
                    and has_amine
                ):
                    print(
                        f"Detected late-stage amide coupling based on functional groups at depth {depth}: {rsmi}"
                    )
                    amide_coupling_found = True
                    return  # Early return once found

        # Only continue traversal if we haven't found an amide coupling yet
        if not amide_coupling_found:
            # Traverse children with incremented depth
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)
                # Early return if amide coupling found in children
                if amide_coupling_found:
                    return

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {amide_coupling_found}")
    return amide_coupling_found
