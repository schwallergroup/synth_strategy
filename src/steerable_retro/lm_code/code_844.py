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
    This function detects a synthetic strategy involving borylation of an aryl bromide
    for subsequent coupling.
    """
    # Track if we found the patterns
    found_borylation = False
    found_bromination = False

    def dfs_traverse(node, depth=0):
        nonlocal found_borylation, found_bromination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for borylation (aryl bromide to boronic acid)
                if not found_borylation:
                    aryl_bromide_pattern = Chem.MolFromSmarts("[c][Br]")
                    boron_reagent_pattern = Chem.MolFromSmarts("[B;$(B(O)(O))]")

                    aryl_bromide_found = False
                    boron_reagent_found = False

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(aryl_bromide_pattern):
                                aryl_bromide_found = True
                            if mol and mol.HasSubstructMatch(boron_reagent_pattern):
                                boron_reagent_found = True
                        except:
                            continue

                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        boronic_acid_pattern = Chem.MolFromSmarts("[c][B;$(B(O)(O))]")
                        if product_mol.HasSubstructMatch(boronic_acid_pattern):
                            if aryl_bromide_found and boron_reagent_found:
                                print(f"Found borylation at depth {depth}")
                                found_borylation = True

                # Check for bromination
                if not found_bromination:
                    brominating_agent_pattern = Chem.MolFromSmarts("[Br]")
                    aryl_pattern = Chem.MolFromSmarts("[c]")

                    brominating_agent_found = False
                    aryl_found = False

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(brominating_agent_pattern):
                                brominating_agent_found = True
                            if mol and mol.HasSubstructMatch(aryl_pattern):
                                aryl_found = True
                        except:
                            continue

                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        aryl_bromide_pattern = Chem.MolFromSmarts("[c][Br]")
                        if product_mol.HasSubstructMatch(aryl_bromide_pattern):
                            if brominating_agent_found and aryl_found:
                                print(f"Found bromination at depth {depth}")
                                found_bromination = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_borylation and found_bromination
