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
    Detects a synthesis strategy involving sequential N-alkylation reactions,
    particularly where multiple nitrogen-containing groups are alkylated in sequence.
    """
    n_alkylation_reactions = 0
    n_alkylation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_reactions, n_alkylation_depths

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Look for patterns indicating N-alkylation
                # This is a simplified check - in practice, you would need more sophisticated analysis
                try:
                    # Check if any reactant contains a halogen (likely alkylating agent)
                    has_halogen = any(re.search(r"Br|Cl|I", r) for r in reactants)

                    # Check if any reactant contains nitrogen
                    has_nitrogen = any("N" in r for r in reactants)

                    # Check if product contains nitrogen
                    product_has_nitrogen = "N" in product

                    if has_halogen and has_nitrogen and product_has_nitrogen:
                        n_alkylation_reactions += 1
                        n_alkylation_depths.append(depth)
                        print(f"Found potential N-alkylation at depth {depth}")
                except:
                    print("Error processing SMILES in reaction")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found sequential N-alkylations
    sequential_alkylation = n_alkylation_reactions >= 2

    # Check if the alkylations are in sequence (adjacent depths)
    if sequential_alkylation and len(n_alkylation_depths) >= 2:
        n_alkylation_depths.sort()
        for i in range(len(n_alkylation_depths) - 1):
            if n_alkylation_depths[i + 1] - n_alkylation_depths[i] == 1:
                print("Found sequential N-alkylations at adjacent depths")
                return True

    print(f"Sequential N-alkylation strategy detected: {sequential_alkylation}")
    print(f"N-alkylation reactions found: {n_alkylation_reactions}")
    return sequential_alkylation
