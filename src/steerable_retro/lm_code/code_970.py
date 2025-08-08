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
    This function detects if the synthetic route includes a halide exchange
    (specifically amine to chloride conversion).
    """
    halide_exchange_found = False

    def dfs_traverse(node):
        nonlocal halide_exchange_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[c][NH2]")
            # Check for chloride in product at the same position
            chloride_pattern = Chem.MolFromSmarts("[c][Cl]")

            amine_present = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(amine_pattern):
                        amine_present = True
                        break
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                chloride_present = product_mol and product_mol.HasSubstructMatch(chloride_pattern)
            except:
                chloride_present = False

            if amine_present and chloride_present:
                print("Amine to chloride conversion detected")
                halide_exchange_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return halide_exchange_found
