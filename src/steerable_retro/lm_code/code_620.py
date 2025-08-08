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
    This function detects if the synthetic route includes a nitro reduction to amine.
    """
    nitro_to_amine_found = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains a nitro group
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
            amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

            reactant_has_nitro = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(nitro_pattern):
                        reactant_has_nitro = True
                        break
                except:
                    continue

            # Check if product has an amine group
            product_has_amine = False
            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(amine_pattern):
                    product_has_amine = True
            except:
                pass

            # If reactant has nitro and product has amine, it's likely a nitro reduction
            if reactant_has_nitro and product_has_amine:
                print("Detected nitro reduction to amine")
                nitro_to_amine_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return nitro_to_amine_found
