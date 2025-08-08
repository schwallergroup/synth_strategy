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
    This function detects if a fluorophenyl group is preserved throughout the synthesis.

    It checks if the fluorophenyl group is present in all non-starting materials in the
    main synthetic pathway. Starting materials (in_stock=True) are excluded from the check.
    """
    # Track if all relevant molecules have fluorophenyl
    all_mols_have_fluorophenyl = True

    # Track if we found any non-starting materials
    found_non_starting_materials = False

    def has_fluorophenyl(mol_smiles):
        """Helper function to check if a molecule has a fluorophenyl group"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Check for fluorine atoms attached to aromatic carbons
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "F":
                # Check if this fluorine is attached to an aromatic carbon
                neighbors = atom.GetNeighbors()
                if neighbors and neighbors[0].GetIsAromatic() and neighbors[0].GetSymbol() == "C":
                    return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal all_mols_have_fluorophenyl, found_non_starting_materials

        if node["type"] == "mol":
            # Skip starting materials (in_stock=True)
            if not node.get("in_stock", False):
                mol_smiles = node.get("smiles", "")
                if mol_smiles:
                    found_non_starting_materials = True

                    # Check if molecule has fluorophenyl group
                    if not has_fluorophenyl(mol_smiles):
                        print(f"Molecule without fluorophenyl group found: {mol_smiles}")
                        all_mols_have_fluorophenyl = False

        elif node["type"] == "reaction":
            # Check if reaction preserves fluorophenyl group
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has fluorophenyl
                reactant_has_fluorophenyl = any(has_fluorophenyl(r) for r in reactants)

                # Check if product has fluorophenyl
                product_has_fluorophenyl = has_fluorophenyl(product)

                # If reactants had fluorophenyl but product doesn't, flag it
                if reactant_has_fluorophenyl and not product_has_fluorophenyl:
                    print(f"Reaction lost fluorophenyl group: {rsmi}")
                    all_mols_have_fluorophenyl = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If we didn't find any non-starting materials, return False
    if not found_non_starting_materials:
        print("No non-starting materials found in the synthesis route")
        return False

    if all_mols_have_fluorophenyl:
        print("Fluorophenyl group is preserved throughout the synthesis")

    return all_mols_have_fluorophenyl
