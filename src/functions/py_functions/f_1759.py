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
    This function detects if the final step in the synthesis is an amide coupling.
    """
    final_step_is_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_coupling

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # The final reaction step is at depth 1 (one step before the final product)
        if node["type"] == "reaction" and depth == 1:
            if "metadata" in node and "rsmi" in node["metadata"]:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing potential final step reaction: {rsmi}")
                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check for amide coupling reaction types
                    amide_reaction_types = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Schotten-Baumann_amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                    ]

                    for reaction_type in amide_reaction_types:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Found amide coupling reaction in final step: {reaction_type}"
                            )
                            final_step_is_amide_coupling = True
                            return

                    # Fallback to functional group checking if reaction type check fails
                    acyl_sources = [
                        "Carboxylic acid",
                        "Acyl halide",
                        "Ester",
                        "Anhydride",
                    ]
                    amine_sources = ["Primary amine", "Secondary amine"]
                    amide_products = [
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                    ]

                    has_acyl_source = any(
                        any(checker.check_fg(fg, r) for r in reactants)
                        for fg in acyl_sources
                    )
                    has_amine_source = any(
                        any(checker.check_fg(fg, r) for r in reactants)
                        for fg in amine_sources
                    )
                    has_amide_product = any(
                        checker.check_fg(fg, product) for fg in amide_products
                    )

                    print(
                        f"FG analysis - Acyl source: {has_acyl_source}, Amine source: {has_amine_source}, Amide product: {has_amide_product}"
                    )

                    if has_acyl_source and has_amine_source and has_amide_product:
                        print(
                            "Found amide coupling based on functional group analysis in final step"
                        )
                        final_step_is_amide_coupling = True

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {final_step_is_amide_coupling}")
    return final_step_is_amide_coupling
