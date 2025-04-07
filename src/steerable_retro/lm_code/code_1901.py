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
    This function detects alpha-chlorination of a ketone.
    """
    alpha_chlorination_detected = False

    def dfs_traverse(node):
        nonlocal alpha_chlorination_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for alpha-chlorination of ketone
                methyl_ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")
                chloro_ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]-[Cl]")

                # Check if methyl ketone is in reactants and chloro ketone is in product
                methyl_ketone_in_reactants = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(methyl_ketone_pattern):
                        methyl_ketone_in_reactants = True
                        break

                product_mol = Chem.MolFromSmiles(product_smiles)
                if (
                    methyl_ketone_in_reactants
                    and product_mol
                    and product_mol.HasSubstructMatch(chloro_ketone_pattern)
                ):
                    print("Alpha-chlorination of ketone detected")
                    alpha_chlorination_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return alpha_chlorination_detected
