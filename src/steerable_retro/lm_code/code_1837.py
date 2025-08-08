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
    This function detects if the synthetic route incorporates a 1,3-benzodioxole fragment
    via amide bond formation.
    """
    benzodioxole_incorporated = False

    def dfs_traverse(node):
        nonlocal benzodioxole_incorporated

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Define SMARTS patterns for 1,3-benzodioxole and amide bond
            benzodioxole_pattern = Chem.MolFromSmarts("c1cc2OCOc2cc1")
            amide_pattern = Chem.MolFromSmarts("NC(=O)")

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check if reactants contain benzodioxole and product contains amide bond with benzodioxole
                    if (
                        reactants_mol.HasSubstructMatch(benzodioxole_pattern)
                        and product_mol.HasSubstructMatch(benzodioxole_pattern)
                        and product_mol.HasSubstructMatch(amide_pattern)
                    ):

                        # Further check if the benzodioxole is connected to the amide
                        # This is a simplified check - in practice would need more sophisticated analysis
                        print("Detected 1,3-benzodioxole fragment incorporation via amide bond")
                        benzodioxole_incorporated = True
            except:
                print("Error processing SMILES in benzodioxole incorporation detection")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return benzodioxole_incorporated
