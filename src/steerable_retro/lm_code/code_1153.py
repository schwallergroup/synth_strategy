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
    Detects if the synthetic route uses a chloroalkyl chain as a leaving group
    for nucleophilic substitution.
    """
    has_chloroalkyl = False

    def dfs_traverse(node):
        nonlocal has_chloroalkyl

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for chloroalkyl pattern in reactants
            chloroalkyl_pattern = Chem.MolFromSmarts("[Cl][CH2][CH2]")

            try:
                for reactant in reactants:
                    if reactant:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(chloroalkyl_pattern):
                            print("Found chloroalkyl leaving group")
                            has_chloroalkyl = True
                            break
            except:
                pass  # Handle parsing errors gracefully

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_chloroalkyl
