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
    Detects if the synthetic route involves formation and subsequent use of a carbamate intermediate.
    """
    carbamate_formation = False
    carbamate_utilization = False

    def dfs_traverse(node):
        nonlocal carbamate_formation, carbamate_utilization

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product) if product else None

            # Check for carbamate formation
            if product_mol:
                carbamate_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])[#8]")
                if product_mol.HasSubstructMatch(carbamate_pattern):
                    print("Found carbamate formation")
                    carbamate_formation = True

            # Check for carbamate utilization (carbamate in reactants)
            for r in reactants:
                reactant_mol = Chem.MolFromSmiles(r) if r else None
                if reactant_mol:
                    carbamate_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])[#8]")
                    if reactant_mol.HasSubstructMatch(carbamate_pattern):
                        print("Found carbamate utilization")
                        carbamate_utilization = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return carbamate_formation and carbamate_utilization
