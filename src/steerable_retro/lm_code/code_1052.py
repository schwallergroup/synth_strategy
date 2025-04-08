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
    This function detects if the synthetic route involves Cbz (benzyloxycarbonyl)
    protection/deprotection of an amine.
    """
    cbz_pattern = Chem.MolFromSmarts("[#6]-[#8]-C(=O)-[#7]")
    cbz_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal cbz_deprotection_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if Cbz group is in reactant but not in product
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(cbz_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and not product_mol.HasSubstructMatch(cbz_pattern):
                            print("Cbz deprotection detected")
                            cbz_deprotection_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cbz_deprotection_detected
