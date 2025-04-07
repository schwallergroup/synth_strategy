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
    This function detects if the final reaction in the synthesis route
    involves an amide formation.
    """
    amide_formation_reactions = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_reactions

        # If we're at a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # We're only interested in the final reaction (depth 1)
            if depth == 1:
                print(
                    f"Examining reaction at depth {depth}: {node['metadata'].get('rsmi', 'No RSMI')}"
                )

                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check if this is an amide formation reaction using the checker
                amide_formation_reactions_list = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Carboxylic acid to amide conversion",
                ]

                for reaction_type in amide_formation_reactions_list:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide formation reaction: {reaction_type}")
                        amide_formation_reactions = True
                        return

                # If no specific reaction type matched, check for functional group changes
                if not amide_formation_reactions:
                    # Check reactants for required functional groups
                    reactants = reactants_str.split(".")
                    has_carbonyl_source = False
                    has_amine = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            has_carbonyl_source = True
                            print(f"Found carbonyl source in reactant: {reactant}")

                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            has_amine = True
                            print(f"Found amine in reactant: {reactant}")

                    # Check if product has amide group
                    has_amide_product = (
                        checker.check_fg("Primary amide", product_str)
                        or checker.check_fg("Secondary amide", product_str)
                        or checker.check_fg("Tertiary amide", product_str)
                    )

                    if has_amide_product:
                        print(f"Found amide in product: {product_str}")

                    # Confirm amide formation
                    if has_carbonyl_source and has_amine and has_amide_product:
                        print("Confirmed amide formation through functional group analysis")
                        amide_formation_reactions = True
                        return

        # Continue DFS traversal with updated depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return amide_formation_reactions
