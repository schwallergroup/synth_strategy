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
    This function detects if the route incorporates multiple chlorine-containing fragments.
    """
    chlorine_fragments = set()

    def dfs_traverse(node):
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            chlorine_pattern = Chem.MolFromSmarts("[Cl]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(chlorine_pattern):
                    # Use a simplified hash of the molecule as identifier
                    chlorine_fragments.add(Chem.MolToSmiles(reactant_mol, isomericSmiles=False))

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    result = len(chlorine_fragments) >= 2
    print(f"Multiple chlorine-containing fragments: {result}")
    print(f"Number of chlorine-containing fragments: {len(chlorine_fragments)}")
    return result
