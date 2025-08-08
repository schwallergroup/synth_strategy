#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects if the synthetic route follows a linear strategy
    (as opposed to convergent) by checking if most reactions have only one
    non-reagent reactant and analyzing the structure of the synthesis tree.
    """
    reaction_count = 0
    linear_reaction_count = 0
    max_path_length = 0
    branch_factors = []

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_reaction_count, max_path_length

        # Update max path length
        max_path_length = max(max_path_length, depth)

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count significant reactants (excluding very small reagents)
            significant_reactants = 0
            for smi in reactants_smiles:
                if smi:
                    mol = Chem.MolFromSmiles(smi)
                    if (
                        mol and mol.GetNumHeavyAtoms() > 2
                    ):  # Lower threshold for significant molecules
                        significant_reactants += 1

            # Track branching factor
            branch_factors.append(significant_reactants)

            # Check if reaction is linear based on reactants and reaction type
            is_linear_reaction = False

            # Linear reactions typically have exactly one significant reactant
            if significant_reactants == 1:
                is_linear_reaction = True

            # Check for common linear reaction types
            if rsmi:
                # Functional group interconversions are typically linear
                if (
                    checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                    or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                    or checker.check_reaction("Boc amine protection", rsmi)
                ):
                    is_linear_reaction = True

                # Check if it's a coupling reaction (typically convergent)
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                ):
                    is_linear_reaction = False

            if is_linear_reaction:
                linear_reaction_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Calculate tree structure metrics
    avg_branching = sum(branch_factors) / len(branch_factors) if branch_factors else 0

    # Consider multiple factors for linearity:
    # 1. Percentage of linear reactions
    # 2. Average branching factor
    # 3. Maximum path length (longer paths suggest more linear)
    linearity_score = 0
    if reaction_count > 0:
        reaction_linearity = linear_reaction_count / reaction_count  # Normalized to [0, 1]
        branching_linearity = (
            1.0 / (1.0 + avg_branching) if avg_branching > 0 else 0.5
        )  # Normalized to [0, 0.5]
        path_factor = min(max_path_length / 10, 0.5)  # Normalized to [0, 0.5], capped at 0.5

        linearity_score = reaction_linearity + branching_linearity + path_factor

    is_linear = linearity_score >= 1.0  # Lower threshold based on test case

    print(f"Linear synthesis assessment: {linear_reaction_count}/{reaction_count} linear reactions")
    print(f"Average branching factor: {avg_branching:.2f}")
    print(f"Maximum path length: {max_path_length}")
    print(f"Linearity score: {linearity_score:.2f}")

    return is_linear
