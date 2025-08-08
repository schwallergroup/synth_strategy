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
    This function detects if the synthetic route uses benzyl protection/deprotection
    of nitrogen atoms as a key strategy.
    """
    # Track if we find both protection and deprotection
    benzyl_protected_intermediates = []
    deprotection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal deprotection_detected

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-benzyl pattern
            nbenzyl_pattern = Chem.MolFromSmarts("[#7]-[#6]-[c]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                has_nbenzyl_product = product_mol.HasSubstructMatch(nbenzyl_pattern)

                # Check each reactant
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        has_nbenzyl_reactant = reactant_mol.HasSubstructMatch(nbenzyl_pattern)

                        # If reactant has N-benzyl but product doesn't, it's a deprotection
                        if has_nbenzyl_reactant and not has_nbenzyl_product:
                            print(f"Detected benzyl deprotection at depth {depth}")
                            deprotection_detected = True

                        # If product has N-benzyl but reactant doesn't, it's a protection
                        if not has_nbenzyl_reactant and has_nbenzyl_product:
                            print(f"Detected benzyl protection at depth {depth}")
                            benzyl_protected_intermediates.append(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found both protection and deprotection
    strategy_detected = deprotection_detected and len(benzyl_protected_intermediates) > 0
    if strategy_detected:
        print("Benzyl protection/deprotection strategy detected in the route")
    return strategy_detected
