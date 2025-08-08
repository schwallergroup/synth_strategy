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
    This function detects a linear synthesis strategy where a core structure
    (like indole) is maintained throughout the synthesis.
    """
    # Initialize tracking variables
    steps_with_indole = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_indole, total_steps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                total_steps += 1
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                    product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                    if product and all(r is not None for r in reactants):
                        # Check for indole core
                        indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
                        if product.HasSubstructMatch(indole_pattern) and any(
                            r.HasSubstructMatch(indole_pattern) for r in reactants
                        ):
                            print(f"Detected indole core preservation at depth {depth}")
                            steps_with_indole += 1

                        # Check for linear synthesis (no more than 2 reactants)
                        if len(reactants) <= 2:
                            print(
                                f"Linear step detected at depth {depth} with {len(reactants)} reactants"
                            )
                except:
                    print("Error processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if the strategy is present
    # At least 3 steps should preserve the indole core in a linear synthesis
    strategy_present = steps_with_indole >= 3 and steps_with_indole / total_steps >= 0.5

    print(
        f"Strategy detection results: {steps_with_indole} out of {total_steps} steps preserve the indole core"
    )

    return strategy_present
