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
    This function detects a synthetic strategy where an acyl group undergoes
    sequential transformations: trifluoroacetyl → carboxylic acid → amide →
    thioamide → thiazole.
    """
    # Track transformations
    found_trifluoroacetyl = False
    found_carboxylic_acid = False
    found_amide = False
    found_thioamide = False
    found_thiazole = False

    def dfs_traverse(node):
        nonlocal found_trifluoroacetyl, found_carboxylic_acid, found_amide, found_thioamide, found_thiazole

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Check for thiazole (depth 0)
                    thiazole_pattern = Chem.MolFromSmarts("c1nc(*)sc1")
                    if product_mol.HasSubstructMatch(thiazole_pattern):
                        found_thiazole = True
                        print("Found thiazole in product")

                    # Check for thioamide (depth 1)
                    thioamide_pattern = Chem.MolFromSmarts("[*]C(=S)[NH2]")
                    if product_mol.HasSubstructMatch(thioamide_pattern):
                        found_thioamide = True
                        print("Found thioamide in product")

                    # Check for amide (depth 2)
                    amide_pattern = Chem.MolFromSmarts("[*]C(=O)[NH2]")
                    if product_mol.HasSubstructMatch(amide_pattern):
                        found_amide = True
                        print("Found amide in product")

                    # Check for carboxylic acid (depth 3)
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[*]C(=O)[OH]")
                    if product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                        found_carboxylic_acid = True
                        print("Found carboxylic acid in product")

                    # Check for trifluoroacetyl (depth 4)
                    trifluoroacetyl_pattern = Chem.MolFromSmarts("[*]C(=O)C(F)(F)F")
                    if product_mol.HasSubstructMatch(trifluoroacetyl_pattern):
                        found_trifluoroacetyl = True
                        print("Found trifluoroacetyl in product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the complete transformation sequence
    return (
        found_trifluoroacetyl
        and found_carboxylic_acid
        and found_amide
        and found_thioamide
        and found_thiazole
    )
