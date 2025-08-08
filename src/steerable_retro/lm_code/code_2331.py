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
    This function detects biaryl formation via Suzuki coupling.
    Looks for reactions where an aryl halide and arylboronic acid form a biaryl product.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant contains boronic acid and another contains halide
            has_boronic_acid = any("[OB]" in r or "B(O)" in r for r in reactants)
            has_halide = any(re.search(r"Br|I|Cl", r) for r in reactants)

            # Check if product has biaryl structure (two connected aromatic rings)
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and has_boronic_acid and has_halide:
                biaryl_pattern = Chem.MolFromSmarts(
                    "[c]!@[c]"
                )  # Two aromatic carbons connected by single bond
                if product_mol.HasSubstructMatch(biaryl_pattern):
                    print("Detected Suzuki coupling for biaryl formation")
                    suzuki_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
