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
    Detects if the synthesis route involves formation of a sulfonamide from an isocyanate.
    """
    isocyanate_found = False
    sulfonamide_formed = False

    def dfs_traverse(node):
        nonlocal isocyanate_found, sulfonamide_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for isocyanate in reactants
                isocyanate_pattern = Chem.MolFromSmarts("[#8]=[#6]=[#7]")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(isocyanate_pattern):
                            isocyanate_found = True
                    except:
                        continue

                # Check for sulfonamide in product
                sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])-[#8]")
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(sulfonamide_pattern):
                        if isocyanate_found:
                            sulfonamide_formed = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Sulfonamide formation from isocyanate: {sulfonamide_formed}")
    return sulfonamide_formed
