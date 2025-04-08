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
    This function detects if the synthesis route involves triflate activation of a phenol.
    """
    triflate_activation_found = False

    def dfs_traverse(node):
        nonlocal triflate_activation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains triflate but reactants contain phenol
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                triflate_pattern = Chem.MolFromSmarts(
                    "[#6]-[#8]-[#16](=[#8])(=[#8])-[#6]([F])([F])[F]"
                )
                if product_mol.HasSubstructMatch(triflate_pattern):
                    # Check if any reactant has a phenol
                    phenol_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            phenol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                            if reactant_mol.HasSubstructMatch(phenol_pattern):
                                phenol_in_reactants = True

                    # If triflate is in product and phenol in reactants, it's triflate activation
                    if phenol_in_reactants:
                        triflate_activation_found = True
                        print("Triflate activation of phenol detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return triflate_activation_found
