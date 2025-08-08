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
    This function detects a synthetic strategy involving biaryl formation via Suzuki coupling.
    """
    suzuki_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for Suzuki coupling pattern
            boronic_acid_pattern = Chem.MolFromSmarts("B(O)(O)c1ccccc1")
            aryl_halide_pattern = Chem.MolFromSmarts("[Cl,Br,I]c1ccccc1")
            biaryl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")

            reactants_with_boronic_acid = any(
                r and r.HasSubstructMatch(boronic_acid_pattern) for r in reactants
            )
            reactants_with_aryl_halide = any(
                r and r.HasSubstructMatch(aryl_halide_pattern) for r in reactants
            )

            if (
                reactants_with_boronic_acid
                and reactants_with_aryl_halide
                and product
                and product.HasSubstructMatch(biaryl_pattern)
            ):
                print(f"Detected Suzuki coupling at depth {depth}")
                suzuki_coupling_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Biaryl formation via Suzuki coupling detected: {suzuki_coupling_detected}")
    return suzuki_coupling_detected
