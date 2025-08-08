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
    This function detects a strategy involving silyl protection and deprotection of amines.
    Specifically, it looks for N-Si bonds being formed or broken in the synthesis route.
    """
    has_silyl_protection = False

    def dfs_traverse(node):
        nonlocal has_silyl_protection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for silyl protection/deprotection
                silyl_pattern = Chem.MolFromSmarts("[#7]-[Si]")

                # Check reactants for silyl groups
                reactants_have_silyl = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(silyl_pattern):
                            reactants_have_silyl = True
                            break
                    except:
                        continue

                # Check product for silyl groups
                product_has_silyl = False
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(silyl_pattern):
                        product_has_silyl = True
                except:
                    pass

                # If silyl group appears or disappears, it's a protection/deprotection
                if (reactants_have_silyl and not product_has_silyl) or (
                    not reactants_have_silyl and product_has_silyl
                ):
                    has_silyl_protection = True
                    print("Detected silyl protection/deprotection")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_silyl_protection
