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
    Detects if the synthesis includes an ester hydrolysis step.
    """
    hydrolysis_found = False

    def dfs_traverse(node):
        nonlocal hydrolysis_found

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester hydrolysis
                product_mol = Chem.MolFromSmiles(product)
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)-[OH]")
                ethyl_ester_pattern = Chem.MolFromSmarts("[C](=O)-O-[CH2]-[CH3]")

                if (
                    product_mol
                    and carboxylic_acid_pattern
                    and product_mol.HasSubstructMatch(carboxylic_acid_pattern)
                ):
                    # Check if reactants include ethyl ester
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if (
                            reactant_mol
                            and ethyl_ester_pattern
                            and reactant_mol.HasSubstructMatch(ethyl_ester_pattern)
                        ):
                            hydrolysis_found = True
                            print("Detected ester hydrolysis step")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return hydrolysis_found
