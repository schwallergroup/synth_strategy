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
    This function detects if the synthesis involves late-stage amide coupling.
    """
    late_stage_amide_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_coupling_detected

        # Check if this is a reaction node at a late stage (depth <= 1)
        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is an amide coupling reaction
                is_amide_coupling = any(
                    [
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        ),
                        checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi),
                        checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                        ),
                        checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi),
                        checker.check_reaction("Ester with primary amine to amide", rsmi),
                        checker.check_reaction("Ester with secondary amine to amide", rsmi),
                        checker.check_reaction("Acylation of primary amines", rsmi),
                        checker.check_reaction("Acylation of secondary amines", rsmi),
                        checker.check_reaction("Schotten-Baumann_amide", rsmi),
                        checker.check_reaction("Acyl chloride with ammonia to amide", rsmi),
                    ]
                )

                # Check for acid in reactants
                has_acid = False
                has_amine = False
                has_acyl_chloride = False
                has_ester = False

                for reactant_smiles in reactants_smiles:
                    if checker.check_fg("Carboxylic acid", reactant_smiles):
                        has_acid = True
                        print(f"Found carboxylic acid in reactant: {reactant_smiles}")
                    if checker.check_fg("Primary amine", reactant_smiles) or checker.check_fg(
                        "Secondary amine", reactant_smiles
                    ):
                        has_amine = True
                        print(f"Found amine in reactant: {reactant_smiles}")
                    if checker.check_fg("Acyl halide", reactant_smiles):
                        has_acyl_chloride = True
                        print(f"Found acyl halide in reactant: {reactant_smiles}")
                    if checker.check_fg("Ester", reactant_smiles):
                        has_ester = True
                        print(f"Found ester in reactant: {reactant_smiles}")

                # Check for amide in product
                has_amide = any(
                    [
                        checker.check_fg("Primary amide", product_smiles),
                        checker.check_fg("Secondary amide", product_smiles),
                        checker.check_fg("Tertiary amide", product_smiles),
                    ]
                )

                if has_amide:
                    print(f"Found amide in product: {product_smiles}")

                # Determine if this is a late-stage amide coupling
                if is_amide_coupling:
                    print(
                        f"Late-stage amide coupling detected (reaction type match) at depth {depth}: {rsmi}"
                    )
                    late_stage_amide_coupling_detected = True
                elif (has_acid or has_acyl_chloride or has_ester) and has_amine and has_amide:
                    print(
                        f"Late-stage amide coupling detected (functional group match) at depth {depth}: {rsmi}"
                    )
                    late_stage_amide_coupling_detected = True
                else:
                    print(
                        f"Not an amide coupling: acid={has_acid}, acyl={has_acyl_chloride}, ester={has_ester}, amine={has_amine}, amide={has_amide}"
                    )

            except Exception as e:
                print(f"Error in late-stage amide coupling detection: {e}")

        # If this is the final product (molecule at depth 0), check its parent reaction
        elif node["type"] == "mol" and depth == 0 and "children" in node and node["children"]:
            # The final product's parent reactions are its children in retrosynthetic tree
            for child in node["children"]:
                if child["type"] == "reaction":
                    dfs_traverse(child, 0)  # Check this reaction as a late-stage reaction

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_amide_coupling_detected
