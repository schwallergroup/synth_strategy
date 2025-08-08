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
    This function detects a synthetic strategy involving O-alkylation of an aromatic hydroxyl
    group, typically used as a protection/modification strategy.
    """
    found_o_alkylation = False
    found_hydroxyl_aromatic = False

    def dfs_traverse(node):
        nonlocal found_o_alkylation, found_hydroxyl_aromatic

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic hydroxyl pattern
                hydroxyl_aromatic_pattern = Chem.MolFromSmarts("[c][OH]")
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(hydroxyl_aromatic_pattern):
                        found_hydroxyl_aromatic = True
                        print("Found aromatic hydroxyl:", reactant)

                # Check for O-alkylation pattern
                if found_hydroxyl_aromatic:
                    alkylated_pattern = Chem.MolFromSmarts("[c][O][C]")
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(alkylated_pattern):
                        found_o_alkylation = True
                        print("Found O-alkylation")

                # Also check for dealkylation (reverse direction)
                if not found_o_alkylation:
                    alkylated_pattern = Chem.MolFromSmarts("[c][O][C]")
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(alkylated_pattern):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(
                                hydroxyl_aromatic_pattern
                            ):
                                found_o_alkylation = True
                                print("Found O-dealkylation (protection removal)")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    strategy_detected = found_o_alkylation and found_hydroxyl_aromatic

    if strategy_detected:
        print("Detected O-alkylation protection strategy")

    return strategy_detected
