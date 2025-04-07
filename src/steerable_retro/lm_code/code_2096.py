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
    Detects if the final step (depth 0) in the synthesis is an N-acylation reaction.
    """
    found_late_acylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_acylation

        print(f"Traversing node of type {node['type']} at depth {depth}")

        # For molecule nodes, just traverse children
        if node["type"] == "mol":
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)
            return

        # For reaction nodes
        if node["type"] == "reaction":
            try:
                # Check if this is a leaf reaction node (no further children) or has FINAL in ID
                is_final_step = (
                    depth == 0
                    or "FINAL"  # First reaction in traversal
                    in node.get("metadata", {}).get("ID", "")
                    or all(  # ID contains FINAL
                        child.get("in_stock", False)
                        for child in node.get("children", [])
                        if child["type"] == "mol"
                    )  # All children are in stock
                )

                if is_final_step:
                    print(f"Examining potential final reaction step")

                    # Extract reactants and product
                    rsmi = node["metadata"].get("rsmi", "")
                    if not rsmi:
                        print("No reaction SMILES found in metadata")
                        return

                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants_smiles}")
                    print(f"Product: {product_smiles}")

                    # Check if this is an N-acylation reaction
                    is_acylation = (
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        )
                        or checker.check_reaction("Acylation of primary amines", rsmi)
                        or checker.check_reaction("Acylation of secondary amines", rsmi)
                        or checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                        )
                        or checker.check_reaction(
                            "Acyl chloride with secondary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction("Ester with primary amine to amide", rsmi)
                        or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                        or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    )

                    if is_acylation:
                        print("Detected N-acylation reaction")
                        found_late_acylation = True
                    else:
                        # Alternative check: look for amide formation
                        product_has_amide = (
                            checker.check_fg("Primary amide", product_smiles)
                            or checker.check_fg("Secondary amide", product_smiles)
                            or checker.check_fg("Tertiary amide", product_smiles)
                        )

                        reactants_have_amide = any(
                            checker.check_fg("Primary amide", r)
                            or checker.check_fg("Secondary amide", r)
                            or checker.check_fg("Tertiary amide", r)
                            for r in reactants_smiles
                        )

                        # Check for amine in reactants
                        reactants_have_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            or checker.check_fg("Aniline", r)
                            for r in reactants_smiles
                        )

                        # Check for acyl source in reactants
                        reactants_have_acyl_source = any(
                            checker.check_fg("Acyl halide", r)
                            or checker.check_fg("Carboxylic acid", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Anhydride", r)
                            for r in reactants_smiles
                        )

                        if (
                            product_has_amide
                            and not reactants_have_amide
                            and reactants_have_amine
                            and reactants_have_acyl_source
                        ):
                            print("Detected amide formation from amine and acyl source")
                            found_late_acylation = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

            # Continue traversal for reaction nodes
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal from the target molecule
    if route["type"] == "mol":
        print("Starting traversal from target molecule")
        dfs_traverse(route)
    else:
        print("Route does not start with a molecule node")

    print(f"Final result: {found_late_acylation}")
    return found_late_acylation
