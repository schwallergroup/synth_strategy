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
    This function detects if the synthetic route involves transformation of an ester to an amide.
    """
    has_transformation = False

    def dfs_traverse(node):
        nonlocal has_transformation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester in reactants
                ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#6]")
                ester_in_reactants = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(ester_pattern):
                        ester_in_reactants = True
                        break

                # Check for amide in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    amide_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#7]")
                    if ester_in_reactants and product_mol.HasSubstructMatch(amide_pattern):
                        has_transformation = True
                        print("Found ester to amide transformation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return has_transformation
