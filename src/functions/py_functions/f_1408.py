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
    This function detects a protection-deprotection sequence of carboxylic acid via methyl ester.
    Looks for carboxylic acid → methyl ester → carboxylic acid pattern.
    """
    # Track if we've seen the pattern
    found_protection = False
    found_deprotection = False

    # SMARTS patterns
    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[O;H1]")
    methyl_ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

    def dfs_traverse(node):
        nonlocal found_protection, found_deprotection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection (carboxylic acid → ester)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(
                        mol and mol.HasSubstructMatch(carboxylic_acid_pattern)
                        for mol in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(methyl_ester_pattern)
                ):
                    found_protection = True
                    print("Found carboxylic acid protection")

                # Check for deprotection (ester → carboxylic acid)
                if (
                    product_mol
                    and any(
                        mol and mol.HasSubstructMatch(methyl_ester_pattern)
                        for mol in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(carboxylic_acid_pattern)
                ):
                    found_deprotection = True
                    print("Found carboxylic acid deprotection")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return found_protection and found_deprotection
