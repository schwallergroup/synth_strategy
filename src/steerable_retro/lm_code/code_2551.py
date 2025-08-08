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
    This function detects if the synthesis maintains chlorinated aromatic rings throughout.
    """
    has_chlorinated_aromatics = False
    all_steps_have_chlorinated_aromatics = True
    step_count = 0

    def dfs_traverse(node):
        nonlocal has_chlorinated_aromatics, all_steps_have_chlorinated_aromatics, step_count

        if node["type"] == "reaction":
            step_count += 1
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            # Check for chlorinated aromatic pattern
            chloro_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")

            reactant_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(products_smiles)

            if chloro_aromatic_pattern is not None:
                step_has_chloro = False

                if reactant_mol is not None and reactant_mol.HasSubstructMatch(
                    chloro_aromatic_pattern
                ):
                    has_chlorinated_aromatics = True
                    step_has_chloro = True

                if product_mol is not None and product_mol.HasSubstructMatch(
                    chloro_aromatic_pattern
                ):
                    has_chlorinated_aromatics = True
                    step_has_chloro = True

                if not step_has_chloro:
                    all_steps_have_chlorinated_aromatics = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Only return true if we have chlorinated aromatics and they appear in all steps
    result = has_chlorinated_aromatics and all_steps_have_chlorinated_aromatics and step_count > 0
    if result:
        print(f"Found chlorinated aromatics maintained throughout {step_count} steps")
    return result
