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
    Detects a linear synthesis route that preserves an aryl bromide throughout
    the synthesis while building molecular complexity.
    """
    is_linear = True
    has_aryl_bromide = False
    aryl_bromide_preserved = False
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, has_aryl_bromide, aryl_bromide_preserved, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Check if this is a branching point (convergent synthesis)
            reactant_count = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    reactant_count += 1

            if reactant_count > 1:
                is_linear = False

            # Check for aryl bromide in product
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                if product:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        aryl_bromide_pattern = Chem.MolFromSmarts("[Br][c]")
                        if product_mol.HasSubstructMatch(aryl_bromide_pattern):
                            has_aryl_bromide = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if aryl bromide is preserved throughout (present in final product)
    if "smiles" in route:
        final_product_mol = Chem.MolFromSmiles(route["smiles"])
        if final_product_mol:
            aryl_bromide_pattern = Chem.MolFromSmarts("[Br][c]")
            if final_product_mol.HasSubstructMatch(aryl_bromide_pattern):
                aryl_bromide_preserved = True

    result = is_linear and has_aryl_bromide and aryl_bromide_preserved and reaction_count >= 3
    print(f"Linear synthesis with preserved aryl bromide: {result}")
    return result
