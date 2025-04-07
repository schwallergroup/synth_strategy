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
    Detects if the final step (depth 0) is an ester hydrolysis to form a carboxylic acid.
    """
    # Track if we found a terminal ester hydrolysis and at what depth
    found_hydrolysis = False
    min_depth_found = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_hydrolysis, min_depth_found

        print(f"Traversing node of type {node['type']} at depth {depth}")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an ester hydrolysis reaction using multiple reaction types
                is_hydrolysis = (
                    checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                )

                # Manual check for ester hydrolysis pattern
                if not is_hydrolysis:
                    # Check if product contains carboxylic acid
                    if checker.check_fg("Carboxylic acid", product_smiles):
                        print("Found carboxylic acid in product")

                        # Check if any reactant contains ester
                        has_ester_reactant = any(
                            checker.check_fg("Ester", reactant) for reactant in reactants_smiles
                        )

                        if has_ester_reactant:
                            print(
                                "Found ester in reactants and carboxylic acid in product - likely an ester hydrolysis"
                            )
                            is_hydrolysis = True

                if is_hydrolysis:
                    print(f"Reaction type matches ester hydrolysis at depth {depth}")

                    # Check if product contains carboxylic acid
                    if checker.check_fg("Carboxylic acid", product_smiles):
                        print("Found carboxylic acid in product")

                        # Check if any reactant contains ester
                        for reactant_smiles in reactants_smiles:
                            if checker.check_fg("Ester", reactant_smiles):
                                print(f"Found ester in reactant: {reactant_smiles}")

                                # Update if this is the lowest depth (most terminal) we've found
                                if depth < min_depth_found:
                                    min_depth_found = depth
                                    found_hydrolysis = True
                                    print(f"Updated terminal hydrolysis at depth {depth}")
                                break
                else:
                    print(f"Reaction at depth {depth} is not an ester hydrolysis")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the hydrolysis was found at depth 1 (terminal reaction)
    print(f"Minimum depth found: {min_depth_found}")
    return found_hydrolysis and min_depth_found == 1
