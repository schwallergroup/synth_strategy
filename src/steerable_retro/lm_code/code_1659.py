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
    Detects if the synthetic route includes amide bond formation using an acid chloride.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acid chloride and amine in reactants, amide in product
            acid_chloride_found = False
            amine_found = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    acid_chloride_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-Cl")
                    amine_pattern = Chem.MolFromSmarts("[#6]-[#7;H2]")

                    if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                        acid_chloride_found = True
                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        amine_found = True

            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                amide_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6](=[#8])-[#6]")
                if product_mol.HasSubstructMatch(amide_pattern):
                    if acid_chloride_found and amine_found:
                        amide_formation_found = True
                        print("Amide formation from acid chloride detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
