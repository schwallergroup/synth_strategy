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
    This function detects a strategy involving sequential phenol O-alkylation reactions
    (methylation, benzylation, and alkylation with an ester).
    """
    phenol_alkylation_count = 0

    def dfs_traverse(node):
        nonlocal phenol_alkylation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol O-alkylation pattern
                # Look for a pattern where a phenol ([OH][c]) becomes an ether ([O][c])
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol:
                            # Check if reactant has a phenol group
                            if react_mol.HasSubstructMatch(Chem.MolFromSmarts("[OH][c]")):
                                # Check if product has a new ether bond where the phenol was
                                if prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[O][c]")):
                                    # Check for specific types of alkylation
                                    if (
                                        prod_mol.HasSubstructMatch(
                                            Chem.MolFromSmarts("[CH3][O][c]")
                                        )
                                        or prod_mol.HasSubstructMatch(
                                            Chem.MolFromSmarts("[c][CH2][O][c]")
                                        )
                                        or prod_mol.HasSubstructMatch(
                                            Chem.MolFromSmarts("[C](=O)[O][CH2][O][c]")
                                        )
                                    ):
                                        phenol_alkylation_count += 1
                                        print(f"Found phenol O-alkylation: {rsmi}")
                                        break
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found at least 2 phenol O-alkylation reactions
    result = phenol_alkylation_count >= 2
    print(
        f"Sequential phenol functionalization detected: {result} (count: {phenol_alkylation_count})"
    )
    return result
