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
    This function detects if the synthetic route involves a linear assembly of three fragments.
    Counts the number of key bond-forming reactions in a linear sequence.
    """
    bond_forming_reactions = 0

    def dfs_traverse(node):
        nonlocal bond_forming_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a bond-forming reaction (more than one reactant)
            if len(reactants) > 1:
                # Check if it's a C-C or C-N bond formation
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    if prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[C]-[C]")
                    ) or prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[C]-[N]")):
                        bond_forming_reactions += 1
                        print(f"Found bond-forming reaction #{bond_forming_reactions}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if we found exactly 2 bond-forming reactions (which would connect 3 fragments)
    return bond_forming_reactions >= 2
