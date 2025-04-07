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
    This function detects the alkylation of a benzofuran core with a dibromoalkane.
    """
    benzofuran_alkylation = False

    def dfs_traverse(node):
        nonlocal benzofuran_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzofuran in reactants
                benzofuran_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc[c]2")
                dibromo_pattern = Chem.MolFromSmarts("Br[CH2][CH2][CH2][CH2]Br")

                has_benzofuran = False
                has_dibromo = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(benzofuran_pattern):
                                has_benzofuran = True
                            if mol.HasSubstructMatch(dibromo_pattern):
                                has_dibromo = True
                    except:
                        continue

                # Check if product has alkylated benzofuran
                if has_benzofuran and has_dibromo:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(benzofuran_pattern):
                            alkylated_pattern = Chem.MolFromSmarts(
                                "c1ccc2c(c1)oc(c2)[CH2][CH2][CH2][CH2]Br"
                            )
                            if prod_mol.HasSubstructMatch(alkylated_pattern):
                                benzofuran_alkylation = True
                                print("Found benzofuran alkylation with dibromoalkane")
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return benzofuran_alkylation
