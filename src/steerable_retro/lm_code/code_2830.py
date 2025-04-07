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
    This function detects if the synthesis includes an ester hydrolysis step.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester in reactants
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
                    if reactant_mol.HasSubstructMatch(ester_pattern):
                        # Check if product has carboxylic acid
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                            if product_mol.HasSubstructMatch(acid_pattern):
                                has_ester_hydrolysis = True
                                print(
                                    f"Detected ester hydrolysis in reaction {node.get('metadata', {}).get('ID', '')}"
                                )
                                break

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_ester_hydrolysis
