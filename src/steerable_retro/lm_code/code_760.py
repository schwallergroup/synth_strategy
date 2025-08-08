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
    Detects if the synthesis uses a convergent approach with thiazole formation as a key step.
    """
    convergent_thiazole_detected = False

    def dfs_traverse(node):
        nonlocal convergent_thiazole_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a convergent step (multiple reactants)
                if len(reactants) >= 2:
                    # Check if product contains thiazole
                    product_mol = Chem.MolFromSmiles(product)
                    thiazole_pattern = Chem.MolFromSmarts("c1scnc1")

                    if product_mol and product_mol.HasSubstructMatch(thiazole_pattern):
                        # Check if this reaction forms the thiazole ring
                        thiazole_in_reactants = False
                        for reactant in reactants:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(thiazole_pattern):
                                thiazole_in_reactants = True
                                break

                        if not thiazole_in_reactants:
                            print("Detected convergent synthesis with thiazole formation")
                            convergent_thiazole_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return convergent_thiazole_detected
