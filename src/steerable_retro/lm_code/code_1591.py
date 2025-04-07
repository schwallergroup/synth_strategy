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
    Detects if the synthesis route employs C-N bond formation as a key step,
    particularly in late-stage fragment coupling.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth <= 1 and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have at least 2 reactants
            if len(reactants) >= 2:
                # Look for patterns indicating C-N bond formation
                # One reactant should have an amine, the other should have a leaving group
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Look for C-N bonds in product that weren't in reactants
                    cn_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                    product_cn_bonds = len(product_mol.GetSubstructMatches(cn_bond_pattern))

                    # Check reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    reactant_cn_bonds = sum(
                        len(mol.GetSubstructMatches(cn_bond_pattern))
                        for mol in reactant_mols
                        if mol
                    )

                    # If product has more C-N bonds than reactants combined, and one reactant has NH2
                    if product_cn_bonds > reactant_cn_bonds and any(
                        mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                    ):
                        print(f"Found C-N bond formation at depth {depth}")
                        result = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
