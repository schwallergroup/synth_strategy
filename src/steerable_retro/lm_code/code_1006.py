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
    This function detects if the synthetic route involves sulfonamide formation.
    """
    sulfonamide_formed = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain sulfonamide
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                reactants_have_sulfonamide = False
                for r in reactants:
                    if r is not None and r.HasSubstructMatch(
                        Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])-[#6]")
                    ):
                        reactants_have_sulfonamide = True
                        break

                # Check if product contains sulfonamide
                product = Chem.MolFromSmiles(product_smiles)
                product_has_sulfonamide = False
                if product is not None and product.HasSubstructMatch(
                    Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])-[#6]")
                ):
                    product_has_sulfonamide = True

                # If sulfonamide is in product but not in reactants, it was formed
                if product_has_sulfonamide and not reactants_have_sulfonamide:
                    sulfonamide_formed = True
                    print("Sulfonamide formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return sulfonamide_formed
