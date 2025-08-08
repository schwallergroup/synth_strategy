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
    Detects if a heterocycle (like pyrazole) is introduced via cross-coupling in the final step.
    """
    final_coupling_detected = False
    heterocycle_introduced = False

    def dfs_traverse(node):
        nonlocal final_coupling_detected, heterocycle_introduced

        if node["type"] == "reaction" and node.get("depth", 0) == 0:  # Final reaction (depth 0)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant is a heterocycle (pyrazole)
            pyrazole_pattern = Chem.MolFromSmarts("[n]1[c][c][n][c]1")
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(pyrazole_pattern):
                        heterocycle_introduced = True
                        print("Found heterocycle in reactants")
                except:
                    continue

            # Check if product has the heterocycle attached to an aromatic ring
            try:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    # Pattern for heterocycle connected to aromatic carbon
                    coupling_pattern = Chem.MolFromSmarts("[c]-[n]1[c][c][n][c]1")
                    if prod_mol.HasSubstructMatch(coupling_pattern):
                        final_coupling_detected = True
                        print("Found heterocycle attached to aromatic ring in product")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return final_coupling_detected and heterocycle_introduced
