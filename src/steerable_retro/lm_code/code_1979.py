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
    Detects if the synthesis route involves N-alkylation using a benzyl chloride
    """
    has_n_alkylation = False

    def dfs_traverse(node):
        nonlocal has_n_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-alkylation with benzyl chloride
                benzyl_chloride_patt = Chem.MolFromSmarts("[c][CH2][Cl]")
                primary_amine_patt = Chem.MolFromSmarts("[NH2][c,C]")
                secondary_amine_patt = Chem.MolFromSmarts("[NH]([c,C])[CH2][c]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(secondary_amine_patt):
                    has_benzyl_chloride = False
                    has_primary_amine = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(benzyl_chloride_patt):
                                has_benzyl_chloride = True
                            if reactant_mol.HasSubstructMatch(primary_amine_patt):
                                has_primary_amine = True

                    if has_benzyl_chloride and has_primary_amine:
                        has_n_alkylation = True
                        print("Found N-alkylation with benzyl chloride")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_n_alkylation
