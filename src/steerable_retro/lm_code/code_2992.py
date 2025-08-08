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
    This function detects if the synthesis involves ester hydrolysis to carboxylic acid.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains carboxylic acid
                product_mol = Chem.MolFromSmiles(product)
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

                # Check if reactants contain ester
                has_ester = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

                        if reactant_mol.HasSubstructMatch(ester_pattern):
                            has_ester = True

                # If product has carboxylic acid and reactant has ester
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(carboxylic_acid_pattern)
                    and has_ester
                ):
                    has_ester_hydrolysis = True
                    print(f"Detected ester hydrolysis at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_ester_hydrolysis
