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
    This function detects alcohol activation via mesylation for subsequent displacement.
    """
    mesylation_detected = False

    def dfs_traverse(node):
        nonlocal mesylation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol to mesylate conversion
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                mesylate_pattern = Chem.MolFromSmarts(
                    "[#6]-[#8]-[#16](=[#8])(=[#8])-[#6]"
                )

                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    p_mol = Chem.MolFromSmiles(product)

                    if r_mol and p_mol:
                        if r_mol.HasSubstructMatch(
                            alcohol_pattern
                        ) and p_mol.HasSubstructMatch(mesylate_pattern):
                            print(f"Mesylation detected: {rsmi}")
                            mesylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return mesylation_detected
