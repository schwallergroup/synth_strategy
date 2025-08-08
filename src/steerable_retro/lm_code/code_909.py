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
    Detects a synthetic strategy involving oxidation of a primary alcohol to a carboxylic acid.
    """
    oxidation_detected = False

    def dfs_traverse(node):
        nonlocal oxidation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for alcohol oxidation to carboxylic acid
                    alcohol_pattern = Chem.MolFromSmarts("[#6][#8;H1]")
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8;H1]")

                    for reactant in reactants:
                        if reactant:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if r_mol and r_mol.HasSubstructMatch(alcohol_pattern):
                                p_mol = Chem.MolFromSmiles(product)
                                if p_mol and p_mol.HasSubstructMatch(carboxylic_acid_pattern):
                                    # Ensure the carbon connected to OH in alcohol is the same as in acid
                                    # This is a simplification - in reality would need atom mapping
                                    oxidation_detected = True
                                    print("Detected alcohol to carboxylic acid oxidation")

                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return oxidation_detected
