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
    This function detects ether cleavage reactions in the synthesis.
    """
    ether_cleavage_found = False

    def dfs_traverse(node):
        nonlocal ether_cleavage_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                products = rsmi.split(">")[-1].split(".")

                # Check for ether in reactants
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        ether_pattern = Chem.MolFromSmarts("[#6]-[OX2]-[#6]")
                        if mol.HasSubstructMatch(ether_pattern):
                            # Check if products contain alcohol
                            for product in products:
                                prod_mol = Chem.MolFromSmiles(product)
                                if prod_mol:
                                    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
                                    if prod_mol.HasSubstructMatch(alcohol_pattern):
                                        ether_cleavage_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Ether cleavage detected: {ether_cleavage_found}")
    return ether_cleavage_found
