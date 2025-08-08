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
    Detects if the synthesis involves the formation of a tertiary alcohol through
    nucleophilic addition to a ketone.
    """
    # Initialize tracking variable
    tertiary_alcohol_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal tertiary_alcohol_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for tertiary alcohol formation
            ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")
            tertiary_alcohol_pattern = Chem.MolFromSmarts("[#6](-[#8;H1])(-[#6])(-[#6])")

            has_ketone = False
            has_tertiary_alcohol = False

            # Check reactants for ketone
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(ketone_pattern):
                        has_ketone = True
                except:
                    continue

            # Check product for tertiary alcohol
            try:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(tertiary_alcohol_pattern):
                    has_tertiary_alcohol = True
            except:
                pass

            if has_ketone and has_tertiary_alcohol:
                tertiary_alcohol_formation = True
                print(f"Detected tertiary alcohol formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Tertiary alcohol formation strategy detected: {tertiary_alcohol_formation}")
    return tertiary_alcohol_formation
