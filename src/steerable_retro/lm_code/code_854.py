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
    This function detects alkene oxidation to aldehyde in the synthetic route.
    """
    alkene_oxidation_detected = False

    def dfs_traverse(node):
        nonlocal alkene_oxidation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check for alkene in reactants
                alkene_pattern = Chem.MolFromSmarts("C=C")
                # Check for aldehyde in product
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                # Check for oxidizing agent
                oxidizing_agent = any("O=" in r for r in reactants)

                alkene_in_reactants = any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(alkene_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                aldehyde_in_product = Chem.MolFromSmiles(product).HasSubstructMatch(
                    aldehyde_pattern
                )

                if alkene_in_reactants and aldehyde_in_product and oxidizing_agent:
                    print("Detected alkene oxidation to aldehyde")
                    alkene_oxidation_detected = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return alkene_oxidation_detected
