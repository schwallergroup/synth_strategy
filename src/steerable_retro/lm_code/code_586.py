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


def main(route):
    """
    This function detects if the synthesis follows a linear strategy rather than convergent.
    It checks if most reactions have only one main reactant contributing to the product scaffold.
    """
    reaction_count = 0
    linear_reaction_count = 0
    max_depth = 0
    branch_count = 0

    def dfs_traverse(node, depth=0, path=None):
        nonlocal reaction_count, linear_reaction_count, max_depth, branch_count

        if path is None:
            path = []

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Count significant reactants (those with more than 5 heavy atoms)
                significant_reactants = 0
                main_scaffold_found = False

                for reactant in reactants_smiles:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol is None:
                        continue

                    # Check if this is a significant reactant
                    if mol.GetNumHeavyAtoms() > 5:
                        significant_reactants += 1

                        # Check if this reactant contains most of the product scaffold
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if (
                            product_mol
                            and mol.GetNumHeavyAtoms() > 0.7 * product_mol.GetNumHeavyAtoms()
                        ):
                            main_scaffold_found = True

                # If only one significant reactant or we found a main scaffold contributor, consider it a linear step
                if significant_reactants <= 1 or main_scaffold_found:
                    linear_reaction_count += 1
                else:
                    # This is a convergent step - check if it's a meaningful branch point
                    if len(node.get("children", [])) > 1:
                        branch_count += 1

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        children = node.get("children", [])
        if len(children) > 0:
            for i, child in enumerate(children):
                # Create a new path for each branch
                new_path = path + [i]
                dfs_traverse(child, depth + 1, new_path)

    dfs_traverse(route)

    # Calculate linearity score based on multiple factors
    linearity_score = 0

    # Factor 1: Percentage of linear reactions
    if reaction_count > 0:
        reaction_linearity = linear_reaction_count / reaction_count
        linearity_score += reaction_linearity * 0.7  # 70% weight

    # Factor 2: Route structure (fewer branches = more linear)
    if reaction_count > 0:
        branch_ratio = 1.0 - (branch_count / reaction_count)
        linearity_score += branch_ratio * 0.3  # 30% weight

    # If at least 70% linearity score, consider it a linear synthesis
    is_linear = linearity_score >= 0.7

    print(
        f"Linear synthesis analysis: {linear_reaction_count}/{reaction_count} reactions are linear"
    )
    print(f"Branch count: {branch_count}, Max depth: {max_depth}")
    print(f"Linearity score: {linearity_score:.2f}")

    return is_linear
