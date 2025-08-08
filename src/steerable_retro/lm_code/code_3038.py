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
    This function detects a strategy involving lactam ring formation.
    It looks for reactions where a lactam ring is formed.
    """
    lactam_formation_found = False

    def dfs_traverse(node):
        nonlocal lactam_formation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains a lactam
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                # SMARTS pattern for lactam: nitrogen connected to carbonyl in a ring
                lactam_pattern = Chem.MolFromSmarts("[#7;R]-[#6;R](=[#8])")
                if product_mol.HasSubstructMatch(lactam_pattern):
                    # Check if reactants don't have the lactam
                    reactants_have_lactam = False
                    for r_smi in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if r_mol and r_mol.HasSubstructMatch(lactam_pattern):
                            reactants_have_lactam = True
                            break

                    if not reactants_have_lactam:
                        print("Lactam formation detected in reaction:", rsmi)
                        lactam_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return lactam_formation_found
