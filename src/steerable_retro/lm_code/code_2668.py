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
    Detects if the synthesis route involves early-stage thiazole ring formation.
    """
    thiazole_formed = False
    thiazole_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal thiazole_formed, thiazole_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains thiazole but reactants don't
            product_mol = Chem.MolFromSmiles(product)
            thiazole_pattern = Chem.MolFromSmarts("c1nc(c)sc1")

            if product_mol and thiazole_pattern:
                product_has_thiazole = product_mol.HasSubstructMatch(thiazole_pattern)

                reactants_have_thiazole = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(thiazole_pattern):
                        reactants_have_thiazole = True
                        break

                if product_has_thiazole and not reactants_have_thiazole:
                    thiazole_formed = True
                    thiazole_depth = depth
                    print(f"Thiazole formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Early stage is defined as depth >= 3
    result = thiazole_formed and thiazole_depth >= 3
    print(f"Early thiazole formation strategy detected: {result}")
    return result
