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
    Detects a synthetic strategy involving early acylation of an aromatic heterocycle.
    """
    early_acylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal early_acylation_found

        if node["type"] == "reaction" and depth >= 3:  # Early steps (depth >= 3)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant is an acylating agent
            acylating_agent_found = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    acyl_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")
                    if mol.HasSubstructMatch(acyl_chloride_pattern):
                        acylating_agent_found = True
                        print(f"Found acylating agent at depth {depth}: {reactant}")

            # Check if product has a new carbonyl group attached to aromatic ring
            if acylating_agent_found:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    aromatic_carbonyl_pattern = Chem.MolFromSmarts("c-[C](=O)")
                    if prod_mol.HasSubstructMatch(aromatic_carbonyl_pattern):
                        early_acylation_found = True
                        print(
                            f"Found early acylation of aromatic system at depth {depth}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return early_acylation_found
