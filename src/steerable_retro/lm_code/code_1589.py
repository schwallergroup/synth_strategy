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
    This function detects if the synthesis uses halogenated aromatics as coupling partners.
    """
    has_halogenated_coupling = False

    def dfs_traverse(node):
        nonlocal has_halogenated_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for multiple halogenated aromatics in a single reaction
            halogenated_aromatics = 0
            for reactant in reactants:
                if reactant:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        if (
                            r_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Br]"))
                            or r_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][I]"))
                            or r_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Cl]"))
                        ):
                            halogenated_aromatics += 1

            if halogenated_aromatics >= 2:
                print("Multiple halogenated aromatics detected in a single reaction")
                has_halogenated_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_halogenated_coupling
