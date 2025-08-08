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
    Detects if the synthesis route includes a thioether formation step (C-S-C bond formation).
    """
    found_thioether_formation = False

    def dfs_traverse(node):
        nonlocal found_thioether_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for thiol in reactants
            thiol_pattern = Chem.MolFromSmarts("[#16;H1]")
            thiol_present = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(thiol_pattern):
                    thiol_present = True
                    break

            # Check for thioether in product
            thioether_pattern = Chem.MolFromSmarts("[#6][#16][#6]")
            product_mol = Chem.MolFromSmiles(product_smiles)
            product_has_thioether = product_mol and product_mol.HasSubstructMatch(thioether_pattern)

            if thiol_present and product_has_thioether:
                found_thioether_formation = True
                print("Found thioether formation step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_thioether_formation
