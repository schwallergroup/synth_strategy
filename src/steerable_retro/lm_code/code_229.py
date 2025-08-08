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
    This function detects benzothiazole synthesis via thiourea cyclization.
    It looks for a thiourea intermediate that cyclizes to form a benzothiazole.
    """
    thiourea_found = False
    benzothiazole_formation = False

    def dfs_traverse(node):
        nonlocal thiourea_found, benzothiazole_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thiourea pattern in reactants
                for reactant in reactants:
                    if reactant:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            thiourea_pattern = Chem.MolFromSmarts("NC(=S)N")
                            if mol.HasSubstructMatch(thiourea_pattern):
                                thiourea_found = True
                                print("Thiourea intermediate found")

                # Check for benzothiazole formation
                if product:
                    mol = Chem.MolFromSmiles(product)
                    if mol:
                        benzothiazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2s1")
                        if mol.HasSubstructMatch(benzothiazole_pattern):
                            # Check if this is a ring formation by comparing with reactants
                            benzothiazole_in_reactants = False
                            for reactant in reactants:
                                if reactant:
                                    r_mol = Chem.MolFromSmiles(reactant)
                                    if r_mol and r_mol.HasSubstructMatch(benzothiazole_pattern):
                                        benzothiazole_in_reactants = True

                            if not benzothiazole_in_reactants:
                                benzothiazole_formation = True
                                print("Benzothiazole formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return thiourea_found and benzothiazole_formation
