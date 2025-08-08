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
    Detects if the synthetic route includes an SNAr reaction with morpholine
    as a nucleophile, typically on an electron-deficient aromatic ring.
    """
    snar_found = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl fluoride or other halide (especially with electron-withdrawing groups)
                aryl_halide_pattern = Chem.MolFromSmarts("[c][F,Cl,Br,I]")
                # Check for morpholine
                morpholine_pattern = Chem.MolFromSmarts("[NH]1CCOC[C]1")
                # Check for aryl-morpholine in product
                aryl_morpholine_pattern = Chem.MolFromSmarts("[c][N]1CCOC[C]1")

                aryl_halide_present = False
                morpholine_present = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            aryl_halide_present = True
                        if mol.HasSubstructMatch(morpholine_pattern):
                            morpholine_present = True

                product_mol = Chem.MolFromSmiles(product)
                aryl_morpholine_present = False
                if product_mol:
                    aryl_morpholine_present = product_mol.HasSubstructMatch(aryl_morpholine_pattern)

                if aryl_halide_present and morpholine_present and aryl_morpholine_present:
                    print(f"Found SNAr with morpholine at depth {depth}")
                    snar_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return snar_found
