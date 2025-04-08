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
    This function detects if the synthetic route involves Weinreb amide formation and subsequent reaction.
    """
    weinreb_formation = False
    weinreb_reaction = False

    def dfs_traverse(node):
        nonlocal weinreb_formation, weinreb_reaction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for Weinreb amide formation
                product = Chem.MolFromSmiles(product_smiles)
                if product:
                    weinreb_pattern = Chem.MolFromSmarts("[C](=[O])[N]([CH3])[O][CH3]")
                    if product.HasSubstructMatch(weinreb_pattern):
                        weinreb_formation = True
                        print(f"Detected Weinreb amide formation: {rsmi}")

                # Check for Weinreb amide reaction
                reactants = Chem.MolFromSmiles(reactants_smiles)
                if reactants:
                    weinreb_pattern = Chem.MolFromSmarts("[C](=[O])[N]([CH3])[O][CH3]")
                    if reactants.HasSubstructMatch(weinreb_pattern):
                        # Check if product is a ketone
                        product = Chem.MolFromSmiles(product_smiles)
                        if product:
                            ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C]")
                            if product.HasSubstructMatch(ketone_pattern):
                                weinreb_reaction = True
                                print(f"Detected Weinreb amide reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Weinreb amide strategy detected: {weinreb_formation and weinreb_reaction}")
    return weinreb_formation and weinreb_reaction
