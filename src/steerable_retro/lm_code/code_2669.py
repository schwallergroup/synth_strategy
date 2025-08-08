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
    Detects if the synthesis route involves late-stage nitrile to hydroxyamidine transformation.
    """
    transformation_detected = False
    transformation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal transformation_detected, transformation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrile in reactants
            nitrile_pattern = Chem.MolFromSmarts("C#N")
            hydroxyamidine_pattern = Chem.MolFromSmarts("[NH]=[C]([NH][OH])")

            nitrile_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if (
                    reactant_mol
                    and nitrile_pattern
                    and reactant_mol.HasSubstructMatch(nitrile_pattern)
                ):
                    nitrile_in_reactants = True
                    break

            # Check for hydroxyamidine in product
            product_mol = Chem.MolFromSmiles(product)
            hydroxyamidine_in_product = False
            if (
                product_mol
                and hydroxyamidine_pattern
                and product_mol.HasSubstructMatch(hydroxyamidine_pattern)
            ):
                hydroxyamidine_in_product = True

            if nitrile_in_reactants and hydroxyamidine_in_product:
                transformation_detected = True
                transformation_depth = depth
                print(f"Nitrile to hydroxyamidine transformation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Late stage is defined as depth <= 1
    result = transformation_detected and transformation_depth <= 1
    print(f"Late nitrile to hydroxyamidine strategy detected: {result}")
    return result
