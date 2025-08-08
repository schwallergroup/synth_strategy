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
    This function detects a synthetic strategy involving alkene formation
    in the synthesis route.
    """
    alkene_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal alkene_formation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for absence of alkene in reactants but presence in product
                for reactant in reactants:
                    if Chem.MolFromSmiles(reactant) is not None:
                        alkene_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
                        reactant_has_alkene = Chem.MolFromSmiles(reactant).HasSubstructMatch(
                            alkene_pattern
                        )

                        if not reactant_has_alkene:
                            if Chem.MolFromSmiles(product) is not None:
                                if Chem.MolFromSmiles(product).HasSubstructMatch(alkene_pattern):
                                    print("Found alkene formation")
                                    alkene_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return alkene_formation_found
