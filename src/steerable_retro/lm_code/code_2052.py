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
    Detects a strategy involving aromatic halogenation (specifically iodination).
    """
    found_aromatic_halogenation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_aromatic_halogenation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic halogenation
                product_mol = Chem.MolFromSmiles(product)
                aryl_iodide_pattern = Chem.MolFromSmarts("c-[I]")

                if product_mol and product_mol.HasSubstructMatch(aryl_iodide_pattern):
                    # Check if reactants don't have the iodide
                    has_iodide_in_reactants = False
                    for reactant in reactants:
                        if "I" in reactant:  # Simple check for iodine
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol and r_mol.HasSubstructMatch(aryl_iodide_pattern):
                                has_iodide_in_reactants = True

                    if not has_iodide_in_reactants:
                        found_aromatic_halogenation = True
                        print(f"Found aromatic iodination at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_aromatic_halogenation
