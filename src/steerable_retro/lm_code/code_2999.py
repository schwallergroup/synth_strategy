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
    Detects a synthesis where the final step is an amide coupling between a carboxylic acid and an amine.
    """
    has_amide_coupling_final_step = False

    def dfs_traverse(node, depth=0):
        nonlocal has_amide_coupling_final_step

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # Check if this is a reaction node at depth 0 (direct child of root)
        # or at depth 1 if the root is a molecule
        is_final_step = False
        if node["type"] == "reaction":
            if depth == 0 or depth == 1:
                is_final_step = True
                print(f"Potential final step found at depth {depth}")

        if node["type"] == "reaction" and is_final_step:
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")
                print(f"Reactants: {reactants_smiles}")
                print(f"Product: {product_smiles}")

                # Check if this is an amide coupling reaction using predefined reaction types
                is_amide_coupling = False

                # List of potential amide coupling reaction types
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected amide coupling reaction: {reaction_type}")
                        is_amide_coupling = True
                        break

                if not is_amide_coupling:
                    # If it's not a known amide coupling reaction, check for the functional groups
                    # Check for carboxylic acid, acyl halide, or ester in reactants
                    has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants_smiles
                    )
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)

                    # Check for amine in reactants (primary, secondary, or aniline)
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants_smiles
                    )

                    # Check for amide in product
                    has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    print(
                        f"Functional group analysis: has_acid={has_acid}, has_acyl_halide={has_acyl_halide}, has_ester={has_ester}, has_amine={has_amine}, has_amide={has_amide}"
                    )

                    # Amide coupling can be from acid+amine, acyl halide+amine, or ester+amine
                    if has_amide and (has_acid or has_acyl_halide or has_ester) and has_amine:
                        print(f"Detected amide coupling based on functional group analysis")
                        is_amide_coupling = True

                if is_amide_coupling:
                    has_amide_coupling_final_step = True
                    print(f"Final step is an amide coupling")
                else:
                    print(f"Not an amide coupling reaction")

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Result: has_amide_coupling_final_step = {has_amide_coupling_final_step}")

    return has_amide_coupling_final_step
