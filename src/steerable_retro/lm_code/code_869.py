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
    Detects if the synthetic route involves multiple protection/deprotection steps,
    indicating a protection group strategy.
    """
    protection_count = 0
    deprotection_count = 0

    def dfs_traverse(node):
        nonlocal protection_count, deprotection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                # Common protection groups
                protection_patterns = [
                    (
                        Chem.MolFromSmarts("[C][Si]([C])([C])[O][C]"),
                        Chem.MolFromSmarts("[O][C]"),
                    ),  # Silyl
                    (Chem.MolFromSmarts("[C](=[O])[O][C]"), Chem.MolFromSmarts("[O][C]")),  # Ester
                    (Chem.MolFromSmarts("[C][O][C]"), Chem.MolFromSmarts("[O][C]")),  # Ether
                    (Chem.MolFromSmarts("[N]([C]=[O])[C]"), Chem.MolFromSmarts("[NH][C]")),  # Amide
                ]

                # Check for protection reactions
                for protected_pattern, unprotected_pattern in protection_patterns:
                    if (
                        product_mol
                        and protected_pattern
                        and product_mol.HasSubstructMatch(protected_pattern)
                    ):
                        # Check if reactants have unprotected functional group
                        for r_mol in reactant_mols:
                            if (
                                r_mol
                                and unprotected_pattern
                                and r_mol.HasSubstructMatch(unprotected_pattern)
                            ):
                                if not r_mol.HasSubstructMatch(protected_pattern):
                                    protection_count += 1
                                    print(f"Found protection reaction: {rsmi}")
                                    break

                # Check for deprotection reactions
                for protected_pattern, unprotected_pattern in protection_patterns:
                    if (
                        product_mol
                        and unprotected_pattern
                        and product_mol.HasSubstructMatch(unprotected_pattern)
                    ):
                        # Check if reactants have protected functional group
                        for r_mol in reactant_mols:
                            if (
                                r_mol
                                and protected_pattern
                                and r_mol.HasSubstructMatch(protected_pattern)
                            ):
                                deprotection_count += 1
                                print(f"Found deprotection reaction: {rsmi}")
                                break
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple protection/deprotection steps are found
    return (protection_count + deprotection_count) >= 2
