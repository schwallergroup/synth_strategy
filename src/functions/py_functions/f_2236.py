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
    Detects if the synthesis uses N-methylation to convert a secondary amine to a tertiary amine
    """
    has_n_methylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_n_methylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-methylation pattern
                has_formaldehyde = False
                has_secondary_amine = False

                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        # Check for formaldehyde or equivalent
                        if reactant.strip() == "C=O" or reactant.strip() == "O=C":
                            has_formaldehyde = True
                        # Check for secondary amine
                        if r_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#7;H1;!$(NC=O)]")
                        ):
                            has_secondary_amine = True

                # Check if product has a tertiary amine with methyl
                if has_secondary_amine:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]-[#7](-[#6])-[#6]")
                    ):
                        print(f"Found N-methylation at depth {depth}")
                        has_n_methylation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_n_methylation
