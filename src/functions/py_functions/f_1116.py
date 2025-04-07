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
    Detects a linear synthesis strategy where a core scaffold (like purine, indole, etc.)
    is preserved while sequential modifications are made to side chains.

    The scaffold may be built in early stages and then preserved in later stages.
    """
    # Track synthesis structure
    reaction_count = 0
    branch_count = 0

    # Common ring scaffolds to check
    common_scaffolds = [
        "purine",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzene",
        "naphthalene",
        "anthracene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
    ]

    # Track scaffold presence at each step
    scaffold_at_step = {}

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, branch_count

        if node["type"] == "reaction":
            reaction_count += 1

            # In retrosynthesis, a reaction with >2 children means convergent synthesis
            if len(node.get("children", [])) > 2:
                branch_count += 1
                print(
                    f"Branch detected at depth {depth} with {len(node.get('children', []))} reactants"
                )

            # Get product and reactants
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for scaffolds in the product
                product_scaffolds = []
                for scaffold in common_scaffolds:
                    if checker.check_ring(scaffold, product_smiles):
                        product_scaffolds.append(scaffold)

                # Store scaffold info for this step
                scaffold_at_step[depth] = {
                    "product": product_smiles,
                    "product_scaffolds": product_scaffolds,
                    "reactants": reactants_smiles,
                    "reactant_scaffolds": [],
                }

                # Check for scaffolds in each reactant
                for reactant in reactants_smiles:
                    reactant_scaffold_list = []
                    for scaffold in common_scaffolds:
                        if checker.check_ring(scaffold, reactant):
                            reactant_scaffold_list.append(scaffold)
                    scaffold_at_step[depth]["reactant_scaffolds"].append(
                        reactant_scaffold_list
                    )

                print(f"Depth {depth}: Product scaffolds: {product_scaffolds}")
                print(
                    f"Depth {depth}: Reactant scaffolds: {scaffold_at_step[depth]['reactant_scaffolds']}"
                )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort depths to analyze from late stage (low depth) to early stage (high depth)
    depths = sorted(scaffold_at_step.keys())

    if len(depths) < 2:
        print("Not enough reaction steps to analyze scaffold preservation")
        return False

    # Find the main scaffold in the final product (lowest depth)
    main_scaffold = None
    has_scaffolds = False

    # Find the first depth where a scaffold appears in the product
    for depth in depths:
        if scaffold_at_step[depth]["product_scaffolds"]:
            main_scaffold = scaffold_at_step[depth]["product_scaffolds"][0]
            has_scaffolds = True
            break

    print(f"Main scaffold identified: {main_scaffold}")

    # Modified scaffold preservation logic
    scaffold_preserved = False
    if main_scaffold:
        # Count consecutive steps where the scaffold appears in products
        consecutive_steps = 0
        current_streak = 0

        for depth in depths:
            if main_scaffold in scaffold_at_step[depth]["product_scaffolds"]:
                current_streak += 1
                consecutive_steps = max(consecutive_steps, current_streak)
            else:
                current_streak = 0

        # Scaffold is preserved if it exists in at least 2 consecutive steps
        scaffold_preserved = consecutive_steps >= 2
        print(f"Scaffold appears in {consecutive_steps} consecutive product steps")

    # A linear synthesis with preserved scaffold should have:
    # 1. No branches (linear)
    # 2. Scaffold preserved in at least 2 consecutive steps
    # 3. At least 2 reactions
    # 4. At least one recognized scaffold
    is_linear = branch_count == 0 and reaction_count > 0

    print(f"Linear synthesis detection:")
    print(f"- Reaction count: {reaction_count}")
    print(f"- Branch count: {branch_count}")
    print(f"- Is linear: {is_linear}")
    print(f"- Has scaffolds: {has_scaffolds}")
    print(f"- Scaffold preserved: {scaffold_preserved}")
    print(f"- Main scaffold: {main_scaffold}")

    return is_linear and scaffold_preserved and reaction_count >= 2 and has_scaffolds
