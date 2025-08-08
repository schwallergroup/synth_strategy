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
    This function detects if the synthetic route involves a Sonogashira coupling (aryl halide + terminal alkyne).
    """
    sonogashira_found = False

    def dfs_traverse(node):
        nonlocal sonogashira_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide and terminal alkyne in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I]")
                terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")

                aryl_halide_found = False
                terminal_alkyne_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_found = True

                    if reactant_mol.HasSubstructMatch(terminal_alkyne_pattern):
                        terminal_alkyne_found = True

                # Check for C-C bond formation in product
                if aryl_halide_found and terminal_alkyne_found:
                    product_mol = Chem.MolFromSmiles(product)
                    aryl_alkyne_pattern = Chem.MolFromSmarts("[c]-[C]#[C]")
                    if product_mol and product_mol.HasSubstructMatch(aryl_alkyne_pattern):
                        print("Sonogashira coupling detected")
                        sonogashira_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return sonogashira_found
