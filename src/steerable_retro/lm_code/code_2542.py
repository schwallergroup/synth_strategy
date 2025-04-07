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
    Detects if the synthesis involves a nitrile to amide transformation.
    """
    nitrile_to_amide_detected = False

    def dfs_traverse(node):
        nonlocal nitrile_to_amide_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check if any reactant contains a nitrile group
                    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check if product has an amide group
                    if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                        for reactant_smiles in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                            if reactant_mol and reactant_mol.HasSubstructMatch(nitrile_pattern):
                                print("Nitrile to amide transformation detected")
                                nitrile_to_amide_detected = True
                                break
                except Exception as e:
                    print(f"Error in nitrile to amide detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitrile_to_amide_detected
