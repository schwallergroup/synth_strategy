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
    Detects if the route involves late-stage introduction of a nitro group.
    Late stage is defined as occurring in the first half of the synthesis (lower depth).
    """
    nitro_introduction_found = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal nitro_introduction_found, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if nitro group is introduced in this reaction
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                # Check if any reactant doesn't have nitro group
                has_nitro_in_reactants = all(
                    Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(nitro_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )

                if not has_nitro_in_reactants:
                    print(f"Nitro group introduction detected at depth {depth}")
                    nitro_introduction_found = True
                    # Check if this is late stage (first half of synthesis)
                    if depth <= max_depth / 2:
                        return True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Determine if nitro introduction was in late stage
    return nitro_introduction_found and max_depth > 0
