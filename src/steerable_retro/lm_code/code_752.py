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
    This function detects if the synthetic route involves ketone → thioamide → amidine conversion.
    """
    ketone_to_thioamide = False
    thioamide_to_amidine = False

    def dfs_traverse(node):
        nonlocal ketone_to_thioamide, thioamide_to_amidine

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ketone to thioamide conversion
            ketone_pattern = Chem.MolFromSmarts("[C](=[O])[#6]")
            thioamide_pattern = Chem.MolFromSmarts("[C](=[S])[N]")

            # Check reactants for ketone
            has_ketone = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(ketone_pattern):
                    has_ketone = True
                    break

            # Check product for thioamide
            product_mol = Chem.MolFromSmiles(product)
            if has_ketone and product_mol and product_mol.HasSubstructMatch(thioamide_pattern):
                ketone_to_thioamide = True
                print(f"Ketone to thioamide conversion detected in reaction: {rsmi}")

            # Check for thioamide to amidine conversion
            amidine_pattern = Chem.MolFromSmarts("[C]([N])=[N]")

            # Check reactants for thioamide
            has_thioamide = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(thioamide_pattern):
                    has_thioamide = True
                    break

            # Check product for amidine
            product_mol = Chem.MolFromSmiles(product)
            if has_thioamide and product_mol and product_mol.HasSubstructMatch(amidine_pattern):
                thioamide_to_amidine = True
                print(f"Thioamide to amidine conversion detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if both conversions are found
    result = ketone_to_thioamide and thioamide_to_amidine
    print(f"Ketone → thioamide → amidine conversion strategy: {result}")
    return result
