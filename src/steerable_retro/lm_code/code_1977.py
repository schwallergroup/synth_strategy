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
    Detects thioether formation via alkylation of thiol.
    """
    thioether_formation_detected = False

    def dfs_traverse(node):
        nonlocal thioether_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thiol pattern
            thiol_pattern = Chem.MolFromSmarts("[SH]")

            # Check for alkyl halide pattern
            alkyl_halide_pattern = Chem.MolFromSmarts("[C]-[Br,Cl,I]")

            # Check for thioether pattern in product
            thioether_pattern = Chem.MolFromSmarts("[#6]-[S]-[#6]")

            has_thiol = False
            has_alkyl_halide = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(thiol_pattern):
                        has_thiol = True
                    if mol.HasSubstructMatch(alkyl_halide_pattern):
                        has_alkyl_halide = True

            product_mol = Chem.MolFromSmiles(product)
            has_thioether = False
            if product_mol and product_mol.HasSubstructMatch(thioether_pattern):
                has_thioether = True

            if has_thiol and has_alkyl_halide and has_thioether:
                print("Thioether formation via thiol alkylation detected")
                thioether_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thioether_formation_detected
