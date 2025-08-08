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
    This function detects if the synthesis includes thiazole ring formation.
    """
    has_thiazole_formation = False

    def dfs_traverse(node):
        nonlocal has_thiazole_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if thiazole is present in the product but not in reactants
                thiazole_pattern = Chem.MolFromSmarts("[#6]1[#16][#6][#6][#7]1")

                # Check reactants
                has_thiazole_reactant = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(thiazole_pattern):
                        has_thiazole_reactant = True
                        break

                # Check product
                prod_mol = Chem.MolFromSmiles(product)
                has_thiazole_product = prod_mol and prod_mol.HasSubstructMatch(thiazole_pattern)

                if not has_thiazole_reactant and has_thiazole_product:
                    has_thiazole_formation = True
                    print("Found thiazole ring formation")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Thiazole formation: {'present' if has_thiazole_formation else 'absent'}")
    return has_thiazole_formation
