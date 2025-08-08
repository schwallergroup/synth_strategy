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
    Detects if the synthesis route includes an aromatic to aliphatic ring conversion,
    specifically phenyl to cyclohexyl transformation.
    """
    found_conversion = False

    def dfs_traverse(node):
        nonlocal found_conversion

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains a phenyl ring
            phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
            cyclohexyl_pattern = Chem.MolFromSmarts("C1CCCCC1")

            reactant_has_phenyl = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(phenyl_pattern):
                        reactant_has_phenyl = True
                        break
                except:
                    continue

            # Check if product contains cyclohexyl
            product_has_cyclohexyl = False
            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(cyclohexyl_pattern):
                    product_has_cyclohexyl = True
            except:
                pass

            if reactant_has_phenyl and product_has_cyclohexyl:
                print("Found aromatic to aliphatic ring conversion")
                found_conversion = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_conversion
