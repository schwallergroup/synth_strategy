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
    This function detects if the synthesis involves a specific sequence of functional group transformations:
    hydroxyl → ether → bromide → nitrogen substitution.
    """
    # Track the transformations observed
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(reactants_mols):
                # Check for hydroxyl → ether transformation
                if (
                    any(mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")) for mol in reactants_mols)
                    and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[O][C]"))
                    and not product_mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))
                ):
                    transformations.append(("hydroxyl_to_ether", depth))
                    print(f"Detected hydroxyl → ether transformation at depth {depth}")

                # Check for ether → bromide transformation (specifically benzylic bromination)
                if any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH3]")) for mol in reactants_mols
                ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2][Br]")):
                    transformations.append(("methyl_to_bromomethyl", depth))
                    print(f"Detected methyl → bromomethyl transformation at depth {depth}")

                # Check for bromide → nitrogen substitution
                if any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2][Br]"))
                    for mol in reactants_mols
                ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c[CH2][n]")):
                    transformations.append(("bromide_to_nitrogen", depth))
                    print(f"Detected bromide → nitrogen substitution at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the transformations follow the expected sequence
    # Sort by depth (higher depth = earlier in synthesis)
    transformations.sort(key=lambda x: x[1], reverse=True)

    # Extract just the transformation types in sequence
    transformation_sequence = [t[0] for t in transformations]

    # Check if our sequence is present (may not be complete)
    expected_sequence = ["hydroxyl_to_ether", "methyl_to_bromomethyl", "bromide_to_nitrogen"]

    # Check if the transformations appear in the correct order
    for i in range(len(expected_sequence) - 1):
        if (
            expected_sequence[i] in transformation_sequence
            and expected_sequence[i + 1] in transformation_sequence
        ):
            if transformation_sequence.index(expected_sequence[i]) > transformation_sequence.index(
                expected_sequence[i + 1]
            ):
                return False

    # Return True if we found at least 2 of the expected transformations in the correct order
    found_count = sum(1 for t in expected_sequence if t in transformation_sequence)
    has_sequence = found_count >= 2

    if has_sequence:
        print(f"Detected functional group transformation sequence: {transformation_sequence}")

    return has_sequence
