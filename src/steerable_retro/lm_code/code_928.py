#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects multiple modifications to a heterocyclic core.
    """
    heterocycle_modifications = 0

    def dfs_traverse(node):
        nonlocal heterocycle_modifications

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for heterocycle in both reactants and products
            reactants_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactants_mol and product_mol:
                # Check for common heterocyclic patterns
                heterocycle_patterns = [
                    "c1ncccc1",  # pyridine
                    "c1ncncc1",  # pyrimidine
                    "c1nccnc1",  # pyrazine
                    "c1ncccn1",  # pyridazine
                    "c1cncnc1",  # pyrimidine (alternative)
                ]

                for pattern in heterocycle_patterns:
                    pattern_mol = Chem.MolFromSmarts(pattern)
                    if reactants_mol.HasSubstructMatch(
                        pattern_mol
                    ) and product_mol.HasSubstructMatch(pattern_mol):
                        # If the same heterocycle is in both reactants and products,
                        # it's likely being modified
                        heterocycle_modifications += 1
                        print(f"Heterocycle modification detected: {pattern}")
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_modifications >= 2  # At least 2 modifications to count as a strategy
