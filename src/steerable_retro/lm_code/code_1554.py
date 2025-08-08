#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if the synthesis involves SNAr reaction with phenol (C-O bond formation where
    carbon is part of aromatic system and was previously attached to a leaving group like Cl).
    """
    has_snar_with_phenol = False

    def dfs_traverse(node):
        nonlocal has_snar_with_phenol

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for phenol pattern in reactants
            phenol_pattern = Chem.MolFromSmarts("c[OH]")
            chloro_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")

            phenol_found = False
            chloro_aromatic_found = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(phenol_pattern):
                        phenol_found = True
                    if mol.HasSubstructMatch(chloro_aromatic_pattern):
                        chloro_aromatic_found = True

            # Check if product has new C-O bond where chlorine was
            if phenol_found and chloro_aromatic_found:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    ether_pattern = Chem.MolFromSmarts("cOc")
                    if product_mol.HasSubstructMatch(ether_pattern):
                        print("Found SNAr reaction with phenol")
                        has_snar_with_phenol = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_snar_with_phenol
