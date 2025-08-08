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
    Detects aromatic bromination as a key functionalization step
    """
    found_aromatic_bromination = False

    def dfs_traverse(node):
        nonlocal found_aromatic_bromination

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic bromination pattern
            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                # SMARTS for aromatic bromide
                aromatic_bromide_pattern = Chem.MolFromSmarts("[c]-[Br]")

                if product_mol.HasSubstructMatch(aromatic_bromide_pattern):
                    # Check if bromine is being added (not already present in all reactants)
                    all_reactants_have_br = True
                    for reactant in reactants:
                        if "Br" not in reactant:
                            all_reactants_have_br = False
                            break

                    if not all_reactants_have_br:
                        found_aromatic_bromination = True
                        print("Found aromatic bromination step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_aromatic_bromination
