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
    Detects a strategy where a heteroatom (S) is exchanged for another (O)
    in a late-stage transformation, particularly in vinyl systems.
    """
    # Initialize tracking variables
    has_heteroatom_exchange = False
    reaction_depths = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal has_heteroatom_exchange, max_depth

        # Track maximum depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Store reaction at its depth
            if depth not in reaction_depths:
                reaction_depths[depth] = []

            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for S to O exchange
                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if (
                    reactant_mol
                    and product_mol
                    and reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C][S]"))
                    and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C][O]"))
                ):
                    has_heteroatom_exchange = True
                    reaction_depths[depth].append("S_to_O_exchange")
                    print(f"Detected S to O exchange at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present and occurs in the late stage (lower depth)
    late_stage_exchange = False
    if has_heteroatom_exchange and max_depth > 0:
        for depth, reactions in reaction_depths.items():
            if "S_to_O_exchange" in reactions and depth <= max_depth // 2:
                late_stage_exchange = True
                print(f"Detected late-stage heteroatom exchange at depth {depth}")

    strategy_present = has_heteroatom_exchange and late_stage_exchange

    if strategy_present:
        print("Detected sequential heteroatom exchange strategy in late stage")
    else:
        print("Sequential heteroatom exchange strategy not detected or not in late stage")

    return strategy_present
