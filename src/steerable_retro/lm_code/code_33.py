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
    This function detects Suzuki coupling reactions forming biaryl systems.
    Looks for reactions where an aryl halide and boronic acid/ester form a biaryl C-C bond.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain aryl halide and boronic acid
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#53,#35]")  # Aryl-I or Aryl-Br
            boronic_acid_pattern = Chem.MolFromSmarts(
                "[c]-[B]([O])[O]"
            )  # Simplified boronic acid pattern

            # Check if product contains biaryl system that wasn't in reactants
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_aryl_halide = False
            has_boronic_acid = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                    if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                if (
                    has_aryl_halide
                    and has_boronic_acid
                    and product_mol
                    and product_mol.HasSubstructMatch(biaryl_pattern)
                ):
                    print("Detected Suzuki coupling for biaryl formation")
                    suzuki_detected = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return suzuki_detected
