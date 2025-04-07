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
    This function detects ester reduction to alcohol in the synthetic route.

    In retrosynthetic analysis:
    - The product (synthetic target) would contain a primary alcohol
    - The reactant (precursor) would contain an ester
    - The reaction would be a reduction of ester to primary alcohol
    """
    print("Starting ester_reduction_alcohol_strategy analysis")
    reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal reduction_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this reaction is an ester reduction using multiple methods

                # Method 1: Check using reaction type
                if checker.check_reaction(
                    "Reduction of ester to primary alcohol", rsmi
                ):
                    print(f"Detected ester reduction reaction by reaction type")
                    reduction_detected = True
                    return

                # Method 2: Check for carboxylic acid reduction
                elif checker.check_reaction(
                    "Reduction of carboxylic acid to primary alcohol", rsmi
                ):
                    print(f"Detected carboxylic acid reduction by reaction type")
                    reduction_detected = True
                    return

                # Method 3: Check for functional group transformation
                ester_in_reactants = any(
                    checker.check_fg("Ester", reactant) for reactant in reactants
                )
                primary_alcohol_in_product = checker.check_fg(
                    "Primary alcohol", product
                )

                if ester_in_reactants and primary_alcohol_in_product:
                    print(f"Detected ester in reactants and primary alcohol in product")

                    # Additional verification: Look for reducing agents or reduction conditions
                    reducing_agents = [
                        "[Al]",
                        "[H]",
                        "LiAlH4",
                        "NaBH4",
                        "BH3",
                        "DIBAL",
                        "LiBH4",
                    ]
                    reagents = rsmi.split(">")[1].split(".")

                    if any(
                        agent in "".join(reagents) for agent in reducing_agents
                    ) or "[HH]" in "".join(reagents):
                        print(f"Confirmed reduction conditions present")
                        reduction_detected = True
                        return

                # Method 4: Look for specific patterns in the transformation
                # Check for ester C(=O)O to CH2OH transformation using atom mapping
                for reactant in reactants:
                    if (
                        "[C" in reactant
                        and "=[O" in reactant
                        and any(
                            x in reactant
                            for x in ["[O:12]", "O[C", "OC", "COC(=O)", "CCO[C]"]
                        )
                    ):
                        if "[CH2" in product and "[OH" in product or "CH2OH" in product:
                            print(f"Detected ester to alcohol pattern in SMILES")
                            reduction_detected = True
                            return

            except Exception as e:
                print(f"Error processing reaction node: {e}")
                print(
                    f"Problematic reaction SMILES: {node.get('metadata', {}).get('rsmi', 'Not available')}"
                )

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Ester reduction to alcohol detected: {reduction_detected}")
    return reduction_detected
