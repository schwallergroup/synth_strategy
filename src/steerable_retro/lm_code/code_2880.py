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
    This function detects a synthetic strategy involving triazole ring formation
    through azide cycloaddition followed by additional heterocycle construction.
    """
    azide_used = False
    triazole_formed = False

    def dfs_traverse(node):
        nonlocal azide_used, triazole_formed

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for azide reactant
                azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=[N-]")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(azide_pattern):
                            azide_used = True
                            print("Azide detected in reactants")
                    except:
                        continue

                # Check for triazole formation
                triazole_pattern = Chem.MolFromSmarts("[n]1[n][n][n]c1")
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(triazole_pattern):
                        # Check if triazole wasn't in reactants
                        triazole_in_reactants = False
                        for reactant in reactants:
                            try:
                                react_mol = Chem.MolFromSmiles(reactant)
                                if react_mol and react_mol.HasSubstructMatch(triazole_pattern):
                                    triazole_in_reactants = True
                            except:
                                continue

                        if not triazole_in_reactants:
                            triazole_formed = True
                            print("Triazole formation detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return azide_used and triazole_formed
