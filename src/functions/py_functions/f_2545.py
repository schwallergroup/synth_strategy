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
    Detects if the synthesis follows a linear pathway with multiple C-N bond formations.
    """
    cn_bond_formations = 0
    linear_synthesis = True

    def dfs_traverse(node):
        nonlocal cn_bond_formations, linear_synthesis

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if synthesis is linear (only one main reactant per step)
                # Expanded list of common reagents to exclude
                main_reactants = [
                    r
                    for r in reactants_smiles
                    if not (
                        r.startswith("O")
                        or r.startswith("[O")
                        or r.startswith("Cl")
                        or r.startswith("[Cl")
                        or r.startswith("Br")
                        or r.startswith("[Br")
                        or r.startswith("I")
                        or r.startswith("[I")
                        or r.startswith("F")
                        or r.startswith("[F")
                        or r.startswith("N")
                        or r.startswith("[N")
                        or r.startswith("S")
                        or r.startswith("[S")
                        or r.startswith("P")
                        or r.startswith("[P")
                        or len(r) < 5  # Short SMILES are likely simple reagents
                    )
                ]

                if len(main_reactants) > 1:
                    print(
                        f"Non-linear step detected: {len(main_reactants)} main reactants: {main_reactants}"
                    )
                    linear_synthesis = False

                try:
                    # Check for C-N bond formation
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol and all(r for r in reactant_mols):
                        # Look for new C-N bonds in the product that weren't in the reactants
                        cn_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]")

                        # Count C-N bonds in reactants
                        reactant_cn_bonds = sum(
                            len(r.GetSubstructMatches(cn_bond_pattern))
                            for r in reactant_mols
                            if r
                        )

                        # Count C-N bonds in product
                        product_cn_bonds = len(
                            product_mol.GetSubstructMatches(cn_bond_pattern)
                        )

                        if product_cn_bonds > reactant_cn_bonds:
                            cn_bond_formations += 1
                            print(
                                f"C-N bond formation detected, total: {cn_bond_formations}"
                            )
                except Exception as e:
                    print(f"Error in C-N bond formation detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Print final status for debugging
    print(
        f"Final status: linear_synthesis={linear_synthesis}, cn_bond_formations={cn_bond_formations}"
    )

    # Return True if synthesis is linear and has at least 2 C-N bond formations
    return linear_synthesis and cn_bond_formations >= 2
