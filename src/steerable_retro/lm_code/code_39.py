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
    This function detects a synthetic strategy involving TMS protection
    of a terminal alkyne.
    """
    terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
    tms_alkyne_pattern = Chem.MolFromSmarts("[C]#[C][Si]([C])([C])[C]")

    has_protection = False

    def dfs_traverse(node):
        nonlocal has_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    products_mol = Chem.MolFromSmiles(products_smiles)

                    if reactants_mol and products_mol:
                        # Check if reactants contain terminal alkyne and products contain TMS-alkyne
                        if reactants_mol.HasSubstructMatch(
                            terminal_alkyne_pattern
                        ) and products_mol.HasSubstructMatch(tms_alkyne_pattern):
                            has_protection = True
                            print(f"Found TMS protection of terminal alkyne: {rsmi}")
                except:
                    print(f"Error processing SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return has_protection
