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
    This function detects a synthetic strategy involving multiple protection/deprotection steps.
    """
    protection_count = 0
    deprotection_count = 0

    def dfs_traverse(node):
        nonlocal protection_count, deprotection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc deprotection
            boc_pattern = Chem.MolFromSmarts("[#6]([#6])([#6])([#6])[#8][#6](=[#8])[#8][#7]")

            # Check for trifluoroacetyl deprotection
            trifluoroacetyl_pattern = Chem.MolFromSmarts("[#6]([#9])([#9])([#9])[#6](=[#8])[#7]")

            # Check for benzyl deprotection
            benzyl_pattern = Chem.MolFromSmarts("c[#6][#7]")

            # Check reactants for protection groups
            reactants_with_protection = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if (
                        mol.HasSubstructMatch(boc_pattern)
                        or mol.HasSubstructMatch(trifluoroacetyl_pattern)
                        or mol.HasSubstructMatch(benzyl_pattern)
                    ):
                        reactants_with_protection += 1

            # Check product for protection groups
            product_mol = Chem.MolFromSmiles(product)
            product_has_protection = False
            if product_mol:
                if (
                    product_mol.HasSubstructMatch(boc_pattern)
                    or product_mol.HasSubstructMatch(trifluoroacetyl_pattern)
                    or product_mol.HasSubstructMatch(benzyl_pattern)
                ):
                    product_has_protection = True

            # If reactants have protection groups but product doesn't, it's a deprotection
            if reactants_with_protection > 0 and not product_has_protection:
                deprotection_count += 1
                print(f"Detected deprotection step: {rsmi}")

            # If product has protection groups but reactants don't, it's a protection
            if product_has_protection and reactants_with_protection == 0:
                protection_count += 1
                print(f"Detected protection step: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if there are multiple protection/deprotection steps
    result = protection_count + deprotection_count >= 2
    print(
        f"Multiple protection/deprotection strategy detected: {result} (Protection: {protection_count}, Deprotection: {deprotection_count})"
    )
    return result
