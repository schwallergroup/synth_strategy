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


def main(route):
    """
    Detects if the synthesis incorporates a fluorinated aromatic group.
    """
    found_fluorinated_aromatic = False

    def dfs_traverse(node, depth=0):
        nonlocal found_fluorinated_aromatic

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Pattern for fluorinated aromatic
                fluoro_aromatic_pattern = Chem.MolFromSmarts("c[F]")

                # Check for fluorinated aromatic in reactants and product
                has_fluoro_aromatic_reactant = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(fluoro_aromatic_pattern):
                        has_fluoro_aromatic_reactant = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                has_fluoro_aromatic_product = (
                    product_mol
                    and product_mol.HasSubstructMatch(fluoro_aromatic_pattern)
                )

                # If fluorinated aromatic is introduced in this step
                if not has_fluoro_aromatic_reactant and has_fluoro_aromatic_product:
                    found_fluorinated_aromatic = True
                    print(
                        f"Detected fluorinated aromatic incorporation at depth {depth}"
                    )
                # Or if it's already present in reactants and preserved in product
                elif has_fluoro_aromatic_reactant and has_fluoro_aromatic_product:
                    found_fluorinated_aromatic = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_fluorinated_aromatic
