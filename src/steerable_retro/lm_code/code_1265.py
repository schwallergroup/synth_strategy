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
    This function detects a strategy involving sequential formation of multiple
    heterocycles (imidazole and pyrazole) in the synthesis route.
    """
    # Track heterocycle formations with their depths
    heterocycle_formations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Convert to RDKit molecules
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(r for r in reactants):
                        # Check for imidazole formation
                        imidazole_pattern = Chem.MolFromSmarts("[nH]1cncc1")
                        reactant_imidazole_count = sum(
                            len(r.GetSubstructMatches(imidazole_pattern)) for r in reactants if r
                        )
                        product_imidazole_count = len(
                            product.GetSubstructMatches(imidazole_pattern)
                        )

                        if product_imidazole_count > reactant_imidazole_count:
                            heterocycle_formations.append(("imidazole", depth))
                            print(f"Detected imidazole formation at depth {depth}")

                        # Check for pyrazole formation
                        pyrazole_pattern = Chem.MolFromSmarts("[nH]1ncc[c]1")
                        reactant_pyrazole_count = sum(
                            len(r.GetSubstructMatches(pyrazole_pattern)) for r in reactants if r
                        )
                        product_pyrazole_count = len(product.GetSubstructMatches(pyrazole_pattern))

                        if product_pyrazole_count > reactant_pyrazole_count:
                            heterocycle_formations.append(("pyrazole", depth))
                            print(f"Detected pyrazole formation at depth {depth}")

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 2 different heterocycle formations
    unique_heterocycles = set(h_type for h_type, _ in heterocycle_formations)

    # The strategy is present if we have at least 2 different heterocycle formations
    strategy_present = len(unique_heterocycles) >= 2

    if strategy_present:
        print(
            f"Detected sequential heterocycle formation strategy with {len(unique_heterocycles)} different heterocycles"
        )
        for h_type, depth in heterocycle_formations:
            print(f"  - {h_type} formation at depth {depth}")

    return strategy_present
