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
    Detects if the route follows a linear strategy of sequential functionalization
    of heterocyclic scaffolds without convergent steps.
    """
    # Track the number of reactants in each step
    multi_reactant_steps = 0
    heterocycle_modifications = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal multi_reactant_steps, heterocycle_modifications, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Count reactants
                if len(reactants) > 1:
                    multi_reactant_steps += 1

                # Check for heterocycle patterns
                heterocycle_patterns = [
                    Chem.MolFromSmarts("c1ccc2c(c1)nccc2"),  # quinoline
                    Chem.MolFromSmarts("N1CCNCC1"),  # piperazine
                    Chem.MolFromSmarts("c1ccncc1"),  # pyridine
                ]

                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol:
                    if any(
                        product_mol.HasSubstructMatch(pattern) for pattern in heterocycle_patterns
                    ):
                        heterocycle_modifications += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Linear strategy has few multi-reactant steps and modifies heterocycles
    is_linear = multi_reactant_steps <= 1 and heterocycle_modifications >= 2

    print(f"Multi-reactant steps: {multi_reactant_steps}")
    print(f"Heterocycle modifications: {heterocycle_modifications}")
    print(f"Maximum depth: {max_depth}")
    print(f"Is linear heterocycle functionalization: {is_linear}")

    return is_linear
