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
    This function detects if the synthetic route involves sequential modifications of an indole scaffold.
    """
    indole_reactions_count = 0

    def dfs_traverse(node):
        nonlocal indole_reactions_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reaction involves indole scaffold
            indole_pattern = Chem.MolFromSmarts(
                "c1ccc2[nX3]ccc2c1"
            )  # Matches both N-H and N-substituted indoles

            if any(
                Chem.MolFromSmiles(r).HasSubstructMatch(indole_pattern) for r in reactants
            ) or Chem.MolFromSmiles(product).HasSubstructMatch(indole_pattern):
                indole_reactions_count += 1
                print(f"Indole scaffold modification detected (total: {indole_reactions_count})")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple reactions involve indole scaffold
    return indole_reactions_count >= 2
