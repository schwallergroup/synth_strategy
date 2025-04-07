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
    Detects the use of multiple different halogen types (F, Cl, Br, I)
    in the synthetic route.
    """
    halogen_types = set()

    def dfs_traverse(node):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for each halogen type
                for halogen, symbol in [
                    ("[F]", "F"),
                    ("[Cl]", "Cl"),
                    ("[Br]", "Br"),
                    ("[I]", "I"),
                ]:
                    pattern = Chem.MolFromSmarts(halogen)
                    if mol.HasSubstructMatch(pattern):
                        halogen_types.add(symbol)

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            for part in rsmi.split(">"):
                for smiles in part.split("."):
                    if smiles:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            for halogen, symbol in [
                                ("[F]", "F"),
                                ("[Cl]", "Cl"),
                                ("[Br]", "Br"),
                                ("[I]", "I"),
                            ]:
                                pattern = Chem.MolFromSmarts(halogen)
                                if mol.HasSubstructMatch(pattern):
                                    halogen_types.add(symbol)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if len(halogen_types) >= 3:
        print(f"Multiple halogen types detected: {', '.join(halogen_types)}")
        return True
    return False
