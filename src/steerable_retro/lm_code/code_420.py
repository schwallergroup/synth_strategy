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
    Detects if the synthesis route involves an SNAr coupling between
    an aniline and a chloropyrimidine.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have two reactants (potential coupling)
            if len(reactants) == 2:
                for r1, r2 in [(0, 1), (1, 0)]:  # Try both combinations
                    reactant1 = Chem.MolFromSmiles(reactants[r1])
                    reactant2 = Chem.MolFromSmiles(reactants[r2])
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant1 and reactant2 and product_mol:
                        # Check for aniline pattern in one reactant
                        aniline_pattern = Chem.MolFromSmarts("[c][NH2]")
                        # Check for chloropyrimidine pattern in the other reactant
                        chloro_pattern = Chem.MolFromSmarts("[c][Cl]")
                        # Check for coupled product pattern
                        coupled_pattern = Chem.MolFromSmarts("[c][NH][c]")

                        if (
                            reactant1.HasSubstructMatch(aniline_pattern)
                            and reactant2.HasSubstructMatch(chloro_pattern)
                            and product_mol.HasSubstructMatch(coupled_pattern)
                        ):
                            print(f"Found SNAr coupling at depth {depth}")
                            result = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return result
