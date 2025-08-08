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
    Detects if the synthetic route uses multiple protection/deprotection steps as a key strategy.
    """
    protection_deprotection_count = 0
    boc_pattern = Chem.MolFromSmarts("[C][C]([C])([C])[O][C](=[O])[N]")
    benzyl_pattern = Chem.MolFromSmarts("c1ccccc1C[O]")

    def dfs_traverse(node):
        nonlocal protection_deprotection_count

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                products_mol = Chem.MolFromSmiles(products_smiles)

                if reactants_mol and products_mol:
                    # Check for Boc deprotection
                    if reactants_mol.HasSubstructMatch(
                        boc_pattern
                    ) and not products_mol.HasSubstructMatch(boc_pattern):
                        protection_deprotection_count += 1
                        print(f"Found Boc deprotection: {rsmi}")

                    # Check for benzyl deprotection
                    if reactants_mol.HasSubstructMatch(
                        benzyl_pattern
                    ) and not products_mol.HasSubstructMatch(benzyl_pattern):
                        protection_deprotection_count += 1
                        print(f"Found benzyl deprotection: {rsmi}")

                    # Check for protection (reverse logic)
                    if not reactants_mol.HasSubstructMatch(
                        boc_pattern
                    ) and products_mol.HasSubstructMatch(boc_pattern):
                        protection_deprotection_count += 1
                        print(f"Found Boc protection: {rsmi}")

                    if not reactants_mol.HasSubstructMatch(
                        benzyl_pattern
                    ) and products_mol.HasSubstructMatch(benzyl_pattern):
                        protection_deprotection_count += 1
                        print(f"Found benzyl protection: {rsmi}")
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total protection/deprotection steps: {protection_deprotection_count}")
    return protection_deprotection_count >= 2
