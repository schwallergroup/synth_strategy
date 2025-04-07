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


def main(route, threshold=2):
    """
    Detects if the synthesis route involves multiple C-N bond formations,
    typically through reductive amination or similar reactions.
    """
    cn_bond_formation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal cn_bond_formation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for C-N bond formation
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Look for patterns suggesting reductive amination or similar C-N bond formation
            # This is a simplified approach - in practice, you'd need more sophisticated reaction classification
            amine_pattern = Chem.MolFromSmarts("[#7]")
            carbonyl_pattern = Chem.MolFromSmarts("[#6]=O")

            has_amine = any(
                mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
            )
            has_carbonyl = any(
                mol and mol.HasSubstructMatch(carbonyl_pattern) for mol in reactant_mols if mol
            )

            # If we have both amine and carbonyl in reactants, and the number of fragments decreases
            # it's likely a C-N bond formation
            if has_amine and has_carbonyl and len(reactants_smiles) > 1:
                print(f"Potential C-N bond formation detected at depth {depth}")
                cn_bond_formation_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Total C-N bond formations detected: {cn_bond_formation_count}")
    return cn_bond_formation_count >= threshold
