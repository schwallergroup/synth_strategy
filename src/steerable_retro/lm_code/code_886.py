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
    This function detects phenol activation via triflation.
    Looks for conversion of phenol to triflate for subsequent coupling.
    """
    triflation_detected = False

    def dfs_traverse(node):
        nonlocal triflation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("c-[OH]")

                # Check for triflate in product
                triflate_pattern = Chem.MolFromSmarts("c-[O][S](=O)(=O)[C]([F])([F])[F]")

                has_phenol = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True
                            break
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    has_triflate = product_mol and product_mol.HasSubstructMatch(triflate_pattern)

                    if has_phenol and has_triflate:
                        print("Detected phenol activation via triflation")
                        triflation_detected = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return triflation_detected
