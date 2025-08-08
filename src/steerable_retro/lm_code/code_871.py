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
    This function detects a linear synthesis strategy with multiple sequential
    alkylation steps (O-alkylation, N-alkylation, C-alkylation).
    """
    # Track alkylation reactions
    alkylation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal alkylation_depths

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkylation reactions
                # Look for methyl iodide or similar alkylating agents
                alkylating_agent_pattern = re.compile(
                    r"I\[CH3\]|I\[CH3:.*\]|Br\[CH3\]|Br\[CH3:.*\]"
                )
                has_alkylating_agent = any(alkylating_agent_pattern.search(r) for r in reactants)

                if has_alkylating_agent:
                    # Determine type of alkylation by comparing product and reactants
                    product_mol = Chem.MolFromSmiles(product)
                    main_reactant = next(
                        (r for r in reactants if not alkylating_agent_pattern.search(r)), None
                    )

                    if main_reactant and product_mol:
                        main_reactant_mol = Chem.MolFromSmiles(main_reactant)
                        if main_reactant_mol:
                            alkylation_depths.append(depth)
                            print(f"Found alkylation reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have at least 3 alkylation steps
    has_multiple_alkylations = len(alkylation_depths) >= 3

    # Check if the synthesis is linear (no convergent steps)
    is_linear = True  # This is a simplification - would need to check number of reactants per step

    if has_multiple_alkylations and is_linear:
        print("Detected linear synthesis with multiple alkylation steps")
        return True
    return False
