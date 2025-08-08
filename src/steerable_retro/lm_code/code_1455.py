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
    This function detects a strategy involving borylation of an aryl halide to prepare for coupling.
    """
    borylation_found = False

    def dfs_traverse(node):
        nonlocal borylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains boronic acid/ester
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    boronic_pattern = Chem.MolFromSmarts("[#6]-[B]([O][#6])[O][#6]")
                    if product_mol.HasSubstructMatch(boronic_pattern):
                        # Check if reactants contain aryl halide
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                aryl_halide_pattern = Chem.MolFromSmarts("c-[Br,I,Cl]")
                                if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                                    borylation_found = True
                                    print("Found borylation of aryl halide")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return borylation_found
