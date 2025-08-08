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
    This function detects furan ring formation in the middle stage of synthesis.
    It looks for reactions where a furan pattern appears in products but not in reactants.
    """
    furan_formation_detected = False
    synthesis_depth = 0
    max_depth = 0

    # First pass to determine synthesis depth
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)
        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)
    mid_synthesis_range = (max_depth // 3, 2 * max_depth // 3)

    def dfs_traverse(node, depth=0):
        nonlocal furan_formation_detected

        if node["type"] == "reaction" and mid_synthesis_range[0] <= depth <= mid_synthesis_range[1]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if furan is formed in this reaction
            furan_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#8]:[#6]1")

            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(furan_pattern):
                # Check if furan was not present in reactants
                furan_in_reactants = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(furan_pattern):
                        furan_in_reactants = True
                        break

                if not furan_in_reactants:
                    print(f"Furan formation detected at depth {depth}")
                    furan_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return furan_formation_detected
