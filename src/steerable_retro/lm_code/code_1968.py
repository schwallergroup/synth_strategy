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
    This function detects if the synthesis involves aromatic bromination.
    """
    bromination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal bromination_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains aromatic bromine
            aromatic_br_pattern = Chem.MolFromSmarts("c[Br]")

            # Check if reactants don't contain aromatic bromine
            has_aromatic_br_in_reactants = False
            for reactant in reactants:
                react_mol = Chem.MolFromSmiles(reactant)
                if react_mol and react_mol.HasSubstructMatch(aromatic_br_pattern):
                    has_aromatic_br_in_reactants = True
                    break

            prod_mol = Chem.MolFromSmiles(product)
            if (
                prod_mol
                and prod_mol.HasSubstructMatch(aromatic_br_pattern)
                and not has_aromatic_br_in_reactants
            ):
                print(f"Aromatic bromination detected at depth {depth}")
                bromination_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return bromination_detected
