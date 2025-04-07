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
    This function detects if the synthesis includes a nitro reduction step.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    if reactant_mol.HasSubstructMatch(nitro_pattern):
                        # Check if product has amine in place of nitro
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            amine_pattern = Chem.MolFromSmarts("[NH2]")
                            if product_mol.HasSubstructMatch(amine_pattern):
                                has_nitro_reduction = True
                                print(
                                    f"Detected nitro reduction in reaction {node.get('metadata', {}).get('ID', '')}"
                                )
                                break

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_nitro_reduction
