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
    This function detects a strategy involving aryl ether formation.
    """
    found_aryl_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_aryl_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains aryl ether but reactants don't
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    aryl_ether_pattern = Chem.MolFromSmarts("c-[O]-[C]")
                    if product_mol.HasSubstructMatch(aryl_ether_pattern):
                        # Check if this is a new formation
                        has_aryl_ether_in_reactants = False
                        for r in reactants:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol and r_mol.HasSubstructMatch(aryl_ether_pattern):
                                has_aryl_ether_in_reactants = True
                                break

                        if not has_aryl_ether_in_reactants:
                            print(f"Found aryl ether formation at depth {depth}")
                            found_aryl_ether_formation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_aryl_ether_formation
