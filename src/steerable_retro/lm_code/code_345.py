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
    This function detects a strategy involving sulfide oxidation to sulfoxide.
    """
    sulfide_oxidation_detected = False

    def dfs_traverse(node):
        nonlocal sulfide_oxidation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfide in reactants
            sulfide_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")

            # Check for sulfoxide in product
            sulfoxide_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])-[#6]")

            try:
                product_mol = Chem.MolFromSmiles(product)

                for reactant in reactants:
                    if reactant:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if (
                            reactant_mol
                            and reactant_mol.HasSubstructMatch(sulfide_pattern)
                            and product_mol
                            and product_mol.HasSubstructMatch(sulfoxide_pattern)
                        ):
                            print("Detected sulfide oxidation in reaction:", rsmi)
                            sulfide_oxidation_detected = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return sulfide_oxidation_detected
