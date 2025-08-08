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
    Detects a strategy involving bromination of a methyl ketone to form a bromoacetyl group.
    """
    methyl_ketone_bromination = False

    def dfs_traverse(node, depth=0):
        nonlocal methyl_ketone_bromination

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for methyl ketone in reactant
            reactant_has_methyl_ketone = False
            for reactant in reactants:
                if (
                    Chem.MolFromSmiles(reactant).GetSubstructMatch(
                        Chem.MolFromSmarts("[#6](=[O])[CH3]")
                    )
                    != ()
                ):
                    reactant_has_methyl_ketone = True

            # Check for bromoacetyl in product
            product_has_bromoacetyl = (
                Chem.MolFromSmiles(product).GetSubstructMatch(
                    Chem.MolFromSmarts("[#6](=[O])[CH2]Br")
                )
                != ()
            )

            if reactant_has_methyl_ketone and product_has_bromoacetyl:
                methyl_ketone_bromination = True
                print(f"Detected methyl ketone bromination at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if methyl_ketone_bromination:
        print("Detected methyl ketone bromination strategy")
        return True

    return False
