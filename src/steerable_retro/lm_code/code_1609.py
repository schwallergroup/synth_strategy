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
    This function detects if the synthesis involves lactam formation via a hydroxamic acid intermediate.
    """
    hydroxamic_acid_found = False
    lactam_from_hydroxamic_acid = False

    def dfs_traverse(node):
        nonlocal hydroxamic_acid_found, lactam_from_hydroxamic_acid

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for hydroxamic acid pattern
            hydroxamic_acid_pattern = Chem.MolFromSmarts("[OH][N]=C")

            # Check for lactam pattern
            lactam_pattern = Chem.MolFromSmarts("[N]C=O")

            product_mol = Chem.MolFromSmiles(product) if product else None

            # Check if product has hydroxamic acid pattern
            if product_mol and product_mol.HasSubstructMatch(hydroxamic_acid_pattern):
                print("Found hydroxamic acid formation")
                hydroxamic_acid_found = True

            # Check if reactant has hydroxamic acid and product has lactam
            if product_mol and product_mol.HasSubstructMatch(lactam_pattern):
                for r in reactants:
                    if r:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(hydroxamic_acid_pattern):
                            print("Found lactam formation from hydroxamic acid")
                            lactam_from_hydroxamic_acid = True
                            break

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return lactam_from_hydroxamic_acid
