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
    This function detects if the synthesis route involves a late-stage amide formation
    (in the final step of the synthesis).
    """
    result = False

    def dfs_traverse(node, depth=0, path=[]):
        nonlocal result

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # Check if this is a reaction node
        if node["type"] == "reaction":
            print(f"Found reaction node at depth: {depth}")

            # Check if this is potentially a final reaction (depth 0 or has no parent reaction)
            if depth == 0 or (len(path) > 0 and path[-1]["type"] == "mol"):
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Checking potential final reaction: {rsmi}")

                    # Check if this is an amide formation reaction using the checker
                    amide_formation_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        "Acyl chloride with ammonia to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with ammonia to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann_amide",
                        "Carboxylic acid to amide conversion",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                    ]

                    for reaction_type in amide_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Late-stage amide formation detected: {reaction_type}")
                            result = True
                            return

                    # If no specific reaction type matched, check for amide formation manually
                    try:
                        reactants_smiles = rsmi.split(">")[0].split(".")
                        product_smiles = rsmi.split(">")[-1]

                        print(f"Reactants: {reactants_smiles}")
                        print(f"Product: {product_smiles}")

                        # Check if product has an amide group
                        product_has_amide = (
                            checker.check_fg("Primary amide", product_smiles)
                            or checker.check_fg("Secondary amide", product_smiles)
                            or checker.check_fg("Tertiary amide", product_smiles)
                        )

                        if product_has_amide:
                            print("Product contains amide group")

                            # Check if reactants have necessary functional groups for amide formation
                            has_acyl_source = False
                            has_amine = False

                            for reactant in reactants_smiles:
                                # Check for acyl sources
                                if (
                                    checker.check_fg("Carboxylic acid", reactant)
                                    or checker.check_fg("Acyl halide", reactant)
                                    or checker.check_fg("Ester", reactant)
                                    or checker.check_fg("Anhydride", reactant)
                                ):
                                    has_acyl_source = True
                                    print(f"Found acyl source: {reactant}")

                                # Check for amine sources
                                if (
                                    checker.check_fg("Primary amine", reactant)
                                    or checker.check_fg("Secondary amine", reactant)
                                    or checker.check_fg("Aniline", reactant)
                                ):
                                    has_amine = True
                                    print(f"Found amine source: {reactant}")

                            # Verify that reactants don't already contain the amide
                            reactants_have_amide = False
                            for reactant in reactants_smiles:
                                if (
                                    checker.check_fg("Primary amide", reactant)
                                    or checker.check_fg("Secondary amide", reactant)
                                    or checker.check_fg("Tertiary amide", reactant)
                                ):
                                    reactants_have_amide = True
                                    print(f"Reactant already contains amide: {reactant}")

                            if has_acyl_source and has_amine and not reactants_have_amide:
                                print("Late-stage amide formation detected in final step")
                                result = True
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")
                else:
                    print("Reaction node missing metadata or rsmi")

        # Add current node to path and traverse children
        new_path = path + [node]
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, new_path)

    # Start traversal from root
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Final result: {result}")
    return result
