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
    Detects if the synthesis involves benzyl protection and deprotection of a phenol.
    """
    has_protection = False
    has_deprotection = False

    def dfs_traverse(node):
        nonlocal has_protection, has_deprotection

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Check for benzyl protection (phenol + benzyl halide -> benzyl ether)
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                benzyl_halide_pattern = Chem.MolFromSmarts("[c][C][Cl,Br,I]")
                benzyl_ether_pattern = Chem.MolFromSmarts("[c][C][O][c]")

                # Check for deprotection (benzyl ether -> phenol)
                if product_mol and product_mol.HasSubstructMatch(phenol_pattern):
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(benzyl_ether_pattern):
                            print("Detected benzyl deprotection")
                            has_deprotection = True

                # Check for protection
                if product_mol and product_mol.HasSubstructMatch(benzyl_ether_pattern):
                    has_phenol = False
                    has_benzyl_halide = False

                    for r_mol in reactant_mols:
                        if r_mol:
                            if r_mol.HasSubstructMatch(phenol_pattern):
                                has_phenol = True
                            if r_mol.HasSubstructMatch(benzyl_halide_pattern):
                                has_benzyl_halide = True

                    if has_phenol and has_benzyl_halide:
                        print("Detected benzyl protection")
                        has_protection = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_protection and has_deprotection
