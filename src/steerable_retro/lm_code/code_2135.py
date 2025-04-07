#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the synthetic route involves reduction of a diester to a diol.
    """
    diester_to_diol_detected = False

    def dfs_traverse(node):
        nonlocal diester_to_diol_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for diester pattern in reactants
                diester_pattern = Chem.MolFromSmarts(
                    "[CX3](=[OX1])[OX2][CX4].[CX3](=[OX1])[OX2][CX4]"
                )
                diol_pattern = Chem.MolFromSmarts("[CX4][OX2H].[CX4][OX2H]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(diester_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(diol_pattern):
                            diester_to_diol_detected = True
                            print("Diester to diol reduction detected")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return diester_to_diol_detected
