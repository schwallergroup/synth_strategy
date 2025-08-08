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
    This function detects if the synthetic route involves ester hydrolysis to form carboxylic acid.
    """
    ester_hydrolysis_detected = False

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for ester in reactants
                reactants = Chem.MolFromSmiles(reactants_smiles)
                if reactants:
                    ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][CH3]")
                    if reactants.HasSubstructMatch(ester_pattern):

                        # Check for carboxylic acid in product
                        product = Chem.MolFromSmiles(product_smiles)
                        if product:
                            acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                            if product.HasSubstructMatch(acid_pattern):
                                ester_hydrolysis_detected = True
                                print(f"Detected ester hydrolysis: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Ester hydrolysis strategy detected: {ester_hydrolysis_detected}")
    return ester_hydrolysis_detected
