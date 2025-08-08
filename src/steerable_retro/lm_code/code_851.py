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
    This function detects Grignard addition to aldehyde forming secondary alcohol.
    """
    grignard_addition_detected = False

    def dfs_traverse(node):
        nonlocal grignard_addition_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check for aldehyde in reactants
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                # Check for Grignard reagent in reactants
                grignard_present = any("Mg" in r for r in reactants)
                # Check for secondary alcohol in product
                sec_alcohol_pattern = Chem.MolFromSmarts("[CH]([OH])[#6]")

                aldehyde_in_reactants = any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(aldehyde_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                sec_alcohol_in_product = product_mol and product_mol.HasSubstructMatch(
                    sec_alcohol_pattern
                )

                if aldehyde_in_reactants and sec_alcohol_in_product and grignard_present:
                    print("Detected Grignard addition to aldehyde")
                    grignard_addition_detected = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return grignard_addition_detected
