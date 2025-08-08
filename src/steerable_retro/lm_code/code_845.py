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
    This function detects a synthetic strategy involving multiple ether formations.
    """
    # Track ether formations
    ether_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal ether_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6][Cl,Br,I]")

                phenol_found = False
                alkyl_halide_found = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(phenol_pattern):
                            phenol_found = True
                        if mol and mol.HasSubstructMatch(alkyl_halide_pattern):
                            alkyl_halide_found = True
                    except:
                        continue

                # Check if product has new ether bond
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    ether_pattern = Chem.MolFromSmarts("[c][O][#6]")
                    if product_mol.HasSubstructMatch(ether_pattern):
                        if phenol_found and alkyl_halide_found:
                            print(f"Found ether formation at depth {depth}")
                            ether_formations += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return ether_formations >= 2
