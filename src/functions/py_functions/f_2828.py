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
    This function detects if the synthesis includes multiple C-N bond formations.
    """
    cn_bond_formations = 0

    def dfs_traverse(node):
        nonlocal cn_bond_formations

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create combined reactants molecule for comparison
            combined_reactants = ".".join(reactants)
            reactants_mol = Chem.MolFromSmiles(combined_reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactants_mol and product_mol:
                # Use reaction difference analysis to detect new C-N bonds
                # This is a simplified approach - in practice, you'd need more sophisticated analysis

                # Check for amide formation
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                if product_mol.HasSubstructMatch(amide_pattern):
                    # Check if this is a new amide bond
                    if not reactants_mol.HasSubstructMatch(amide_pattern):
                        cn_bond_formations += 1
                        print(
                            f"Detected C-N bond formation (amide) in reaction {node.get('metadata', {}).get('ID', '')}"
                        )

                # Check for aromatic C-N bond formation
                if "Cl" in "".join(reactants) and "[NH" in "".join(reactants):
                    # This is a simplification - in practice, you'd need to analyze the actual reaction
                    aromatic_cn_pattern = Chem.MolFromSmarts("c[NH]")
                    if product_mol.HasSubstructMatch(aromatic_cn_pattern):
                        cn_bond_formations += 1
                        print(
                            f"Detected aromatic C-N bond formation in reaction {node.get('metadata', {}).get('ID', '')}"
                        )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cn_bond_formations >= 2  # Return True if at least 2 C-N bond formations
