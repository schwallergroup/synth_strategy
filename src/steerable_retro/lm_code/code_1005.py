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
    This function detects if the synthetic route involves thiazole ring formation.
    """
    thiazole_formed = False

    def dfs_traverse(node):
        nonlocal thiazole_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain thiazole
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                reactants_have_thiazole = False
                for r in reactants:
                    if r is not None and r.HasSubstructMatch(Chem.MolFromSmarts("c1scnc1")):
                        reactants_have_thiazole = True
                        break

                # Check if product contains thiazole
                product = Chem.MolFromSmiles(product_smiles)
                product_has_thiazole = False
                if product is not None and product.HasSubstructMatch(Chem.MolFromSmarts("c1scnc1")):
                    product_has_thiazole = True

                # If thiazole is in product but not in reactants, it was formed
                if product_has_thiazole and not reactants_have_thiazole:
                    thiazole_formed = True
                    print("Thiazole ring formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return thiazole_formed
