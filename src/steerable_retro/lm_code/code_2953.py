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
    This function detects the reduction of a carbonyl group to an amine
    in the synthesis route.
    """
    carbonyl_reduction_detected = False

    def dfs_traverse(node):
        nonlocal carbonyl_reduction_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Carbonyl pattern
                carbonyl_pattern = Chem.MolFromSmarts("[#6](=[O])")

                # Amine pattern
                amine_pattern = Chem.MolFromSmarts("[#6][NH2]")

                try:
                    # Check if any reactant has carbonyl group
                    carbonyl_present = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(carbonyl_pattern):
                            carbonyl_present = True
                            break

                    # Check if product has amine group
                    product_mol = Chem.MolFromSmiles(product)
                    amine_present = product_mol and product_mol.HasSubstructMatch(amine_pattern)

                    if carbonyl_present and amine_present:
                        print("Detected carbonyl reduction to amine")
                        carbonyl_reduction_detected = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return carbonyl_reduction_detected
