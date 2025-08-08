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
    This function detects if the synthetic route involves Boc protection/deprotection steps.
    """
    boc_protection_count = 0
    boc_deprotection_count = 0

    def dfs_traverse(node):
        nonlocal boc_protection_count, boc_deprotection_count

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection
            boc_pattern = Chem.MolFromSmarts("[#6][#8][C](=[O])[N]")

            # Check if product has Boc group but reactants don't
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                has_boc_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                        has_boc_in_reactants = True
                        break

                if not has_boc_in_reactants:
                    boc_protection_count += 1
                    print(f"Boc protection detected in reaction: {rsmi}")

            # Check for Boc deprotection
            reactants_have_boc = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                    reactants_have_boc = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            if (
                reactants_have_boc
                and product_mol
                and not product_mol.HasSubstructMatch(boc_pattern)
            ):
                boc_deprotection_count += 1
                print(f"Boc deprotection detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if at least one Boc protection and one deprotection are found
    result = boc_protection_count > 0 and boc_deprotection_count > 0
    print(
        f"Boc protection/deprotection strategy: {result} (Protection: {boc_protection_count}, Deprotection: {boc_deprotection_count})"
    )
    return result
