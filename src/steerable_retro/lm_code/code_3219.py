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
    Detects if the synthesis route involves thioether (C-S-C) bond formation
    as a key fragment coupling strategy.
    """
    thioether_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal thioether_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thioether formation
                thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
                thiol_pattern = Chem.MolFromSmarts("[#6]-[#16H]")
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8H]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(thioether_pattern):
                    # Check if reactants contain thiol and alcohol
                    has_thiol = False
                    has_alcohol = False

                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            if r_mol.HasSubstructMatch(thiol_pattern):
                                has_thiol = True
                            if r_mol.HasSubstructMatch(alcohol_pattern):
                                has_alcohol = True

                    if has_thiol and has_alcohol:
                        print(f"Thioether formation detected at depth {depth}")
                        thioether_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return thioether_formation_detected
