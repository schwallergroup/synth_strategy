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
    Detects if the synthesis route involves protection and deprotection of phenol groups
    """
    has_protection = False
    has_deprotection = False

    def dfs_traverse(node):
        nonlocal has_protection, has_deprotection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for MOM deprotection (MOM-protected phenol -> phenol)
                mom_phenol_patt = Chem.MolFromSmarts("[c][O][CH2][O][CH3]")
                phenol_patt = Chem.MolFromSmarts("[c][OH]")

                product_mol = Chem.MolFromSmiles(product)

                # Check for deprotection
                if product_mol and product_mol.HasSubstructMatch(phenol_patt):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            mom_phenol_patt
                        ):
                            has_deprotection = True
                            print("Found phenol deprotection step")

                # Check for protection
                if product_mol and product_mol.HasSubstructMatch(mom_phenol_patt):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(phenol_patt):
                            has_protection = True
                            print("Found phenol protection step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_protection and has_deprotection
