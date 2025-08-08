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
    This function detects synthesis strategies centered around a thiophene core.
    """
    thiophene_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal thiophene_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Thiophene pattern
                thiophene_pattern = Chem.MolFromSmarts("[c]1[c][s][c][c]1")

                # Check if thiophene is involved in the reaction
                thiophene_involved = False

                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(thiophene_pattern):
                        for reactant in reactants:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(thiophene_pattern):
                                thiophene_involved = True
                                break

                        if thiophene_involved:
                            thiophene_reactions += 1
                            print(f"Thiophene modification detected at depth {depth}")
                except:
                    pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if multiple thiophene-involving reactions are detected
    return thiophene_reactions >= 2
