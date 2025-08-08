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
    This function detects whether the aromatic core structure is preserved throughout the synthesis.
    """
    final_product_smiles = None
    starting_materials = []

    # First pass: collect final product and starting materials
    def collect_molecules(node):
        nonlocal final_product_smiles, starting_materials

        if node.get("type") == "mol":
            if not node.get("children"):  # Leaf node = starting material
                starting_materials.append(node.get("smiles", ""))
            elif final_product_smiles is None:  # Root node = final product
                final_product_smiles = node.get("smiles", "")

        # Traverse children
        for child in node.get("children", []):
            collect_molecules(child)

    # Start collection
    collect_molecules(route)

    # Check if aromatic core is preserved
    if final_product_smiles:
        final_mol = Chem.MolFromSmiles(final_product_smiles)

        # Define core aromatic patterns to check
        isoxazoline_pattern = Chem.MolFromSmarts("[c]1[c][c][c]([C]2=[N][O][C][C]2)[c][c]1")
        dichlorophenyl_pattern = Chem.MolFromSmarts("[c]1[c]([Cl])[c][c]([Cl])[c][c]1")

        core_preserved = False
        if final_mol:
            # Check if these patterns exist in both starting materials and final product
            final_has_isoxazoline = final_mol.HasSubstructMatch(isoxazoline_pattern)
            final_has_dichlorophenyl = final_mol.HasSubstructMatch(dichlorophenyl_pattern)

            starting_has_isoxazoline = False
            starting_has_dichlorophenyl = False

            for sm in starting_materials:
                sm_mol = Chem.MolFromSmiles(sm)
                if sm_mol:
                    if sm_mol.HasSubstructMatch(isoxazoline_pattern):
                        starting_has_isoxazoline = True
                    if sm_mol.HasSubstructMatch(dichlorophenyl_pattern):
                        starting_has_dichlorophenyl = True

            core_preserved = (final_has_isoxazoline and starting_has_isoxazoline) or (
                final_has_dichlorophenyl and starting_has_dichlorophenyl
            )

    print(f"Aromatic core preservation: {core_preserved}")
    return core_preserved
