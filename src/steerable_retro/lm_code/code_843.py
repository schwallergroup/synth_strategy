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
    This function detects a synthetic strategy involving Suzuki coupling for biaryl formation.
    """
    # Track if we found the pattern
    found_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid and aryl halide in reactants
                boronic_acid_pattern = Chem.MolFromSmarts("[c][B;$(B(O)(O))]")
                aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

                boronic_acid_found = False
                aryl_halide_found = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                            boronic_acid_found = True
                        if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                            aryl_halide_found = True
                    except:
                        continue

                # Check if product has new biaryl bond
                if boronic_acid_found and aryl_halide_found:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                        if product_mol.HasSubstructMatch(biaryl_pattern):
                            print(f"Found Suzuki coupling at depth {depth}")
                            found_suzuki = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_suzuki
