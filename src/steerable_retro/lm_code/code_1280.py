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
    Detects synthesis strategy involving reduction of a nitro group to an amine,
    particularly in the context of preparing an aminophenol for further reactions.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Check if this is a reaction node with metadata
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro pattern in reactants
                nitro_pattern = Chem.MolFromSmarts("[#7](=[#8])[#8]")

                # Check for amine pattern in product
                amine_pattern = Chem.MolFromSmarts("[#7;H2]-[#6]")

                nitro_found = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(nitro_pattern):
                            nitro_found = True
                            print("Found nitro group in reactant")
                    except:
                        continue

                # Check if product has amine group
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(amine_pattern) and nitro_found:
                        found_pattern = True
                        print("Found nitro reduction to amine")
                except:
                    pass

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if found_pattern:
        print("Detected nitro reduction to amine strategy")

    return found_pattern
