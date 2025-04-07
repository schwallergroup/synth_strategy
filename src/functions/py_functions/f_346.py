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
    This function detects a strategy involving ester reduction to primary alcohol.
    """
    ester_reduction_detected = False

    def dfs_traverse(node):
        nonlocal ester_reduction_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester in reactants
            ester_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8]-[#6]")

            # Check for primary alcohol in product
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")

            try:
                product_mol = Chem.MolFromSmiles(product)

                for reactant in reactants:
                    if reactant:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if (
                            reactant_mol
                            and reactant_mol.HasSubstructMatch(ester_pattern)
                            and product_mol
                            and product_mol.HasSubstructMatch(alcohol_pattern)
                        ):
                            print("Detected ester reduction in reaction:", rsmi)
                            ester_reduction_detected = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ester_reduction_detected
