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
    Detects a synthetic strategy involving benzyl alcohol → benzyl chloride → benzyl nitrile transformation sequence.
    """
    # Track if we found the required transformations
    found_alcohol_to_chloride = False
    found_chloride_to_nitrile = False

    def dfs_traverse(node):
        nonlocal found_alcohol_to_chloride, found_chloride_to_nitrile

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for benzyl alcohol to benzyl chloride transformation
            alcohol_pattern = Chem.MolFromSmarts("[c][CH2][OH]")
            chloride_pattern = Chem.MolFromSmarts("[c][CH2][Cl]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(m and m.HasSubstructMatch(alcohol_pattern) for m in reactant_mols)
                and product_mol.HasSubstructMatch(chloride_pattern)
            ):
                found_alcohol_to_chloride = True
                print("Found benzyl alcohol to benzyl chloride transformation")

            # Check for benzyl chloride to benzyl nitrile transformation
            nitrile_pattern = Chem.MolFromSmarts("[c][CH2][C]#[N]")

            if (
                product_mol
                and any(m and m.HasSubstructMatch(chloride_pattern) for m in reactant_mols)
                and product_mol.HasSubstructMatch(nitrile_pattern)
            ):
                found_chloride_to_nitrile = True
                print("Found benzyl chloride to benzyl nitrile transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both transformations were found
    return found_alcohol_to_chloride and found_chloride_to_nitrile
