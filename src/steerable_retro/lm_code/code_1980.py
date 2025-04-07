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
    Detects if the synthesis route involves O-methylation of a phenol
    """
    has_o_methylation = False

    def dfs_traverse(node):
        nonlocal has_o_methylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for O-methylation of phenol
                phenol_patt = Chem.MolFromSmarts("[c][OH]")
                methoxy_patt = Chem.MolFromSmarts("[c][O][CH3]")
                methyl_donor_patt = Chem.MolFromSmarts("[CH3][I,Br]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(methoxy_patt):
                    has_phenol = False
                    has_methyl_donor = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(phenol_patt):
                                has_phenol = True
                            if reactant_mol.HasSubstructMatch(methyl_donor_patt):
                                has_methyl_donor = True

                    if has_phenol and has_methyl_donor:
                        has_o_methylation = True
                        print("Found O-methylation of phenol")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_o_methylation
