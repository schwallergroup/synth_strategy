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
    Detects if the synthesis route involves an ester hydrolysis step.
    """
    found_ester_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, product has ester and reactant has carboxylic acid
                ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
                acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]")

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and product_mol.HasSubstructMatch(ester_pattern):
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if not reactant_mol:
                                continue

                            if reactant_mol.HasSubstructMatch(acid_pattern):
                                print(f"Found ester hydrolysis at depth {depth}")
                                found_ester_hydrolysis = True
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_ester_hydrolysis
