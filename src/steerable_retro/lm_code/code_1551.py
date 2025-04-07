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
    Detects synthesis routes that use late-stage amide coupling as a key strategy.
    """
    found_late_amide_coupling = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_late_amide_coupling, max_depth

        # Track maximum depth to determine what's "late-stage"
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide coupling
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Look for carboxylic acid in reactants
                acid_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#8;H1]")
                amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(NC=O)]")
                amide_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#7]")

                has_acid = any(mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols)
                has_amine = any(mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols)
                has_amide_product = product_mol.HasSubstructMatch(amide_pattern)

                if (
                    has_acid and has_amine and has_amide_product and depth <= 1
                ):  # Depth <= 1 means late-stage
                    found_late_amide_coupling = True
                    print(f"Found late-stage amide coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide coupling strategy detected: {found_late_amide_coupling}")
    return found_late_amide_coupling
