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
    This function detects if the synthesis uses TBDMS protection of an alcohol
    early in the synthesis (higher depth).
    """
    tbdms_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal tbdms_protection_found

        if node["type"] == "reaction" and depth >= 2:  # Early in synthesis (higher depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a protection reaction
                tbdms_pattern = Chem.MolFromSmarts("[O][Si]([C])([C])[C]")
                cl_si_pattern = Chem.MolFromSmarts("Cl[Si]")

                # Check if TBDMS appears in product but not in all reactants
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(tbdms_pattern):
                    # Check if one of the reactants is a chlorosilane
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(cl_si_pattern):
                            tbdms_protection_found = True
                            print(f"TBDMS protection found at depth {depth}")
                            break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return tbdms_protection_found
