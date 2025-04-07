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
    This function detects if the final step in the synthesis is an amide coupling.
    """
    # Track if we found late-stage amide coupling
    found_late_amide = False

    # SMARTS patterns
    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[O;H1]")
    amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")

    def dfs_traverse(node):
        nonlocal found_late_amide

        if node["type"] == "reaction" and node.get("depth", 0) == 0:  # Check if it's the final step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                has_acid = any(
                    mol and mol.HasSubstructMatch(carboxylic_acid_pattern) for mol in reactant_mols
                )
                has_amine = any(
                    mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                )
                has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if has_acid and has_amine and has_amide:
                    found_late_amide = True
                    print("Found late-stage amide coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_late_amide
