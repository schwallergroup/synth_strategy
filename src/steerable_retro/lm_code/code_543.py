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
    Detects a convergent synthesis strategy that joins two complex heterocyclic fragments.
    """
    convergent_step_found = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_step_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a convergent step joining heterocycles
            if len(reactants) >= 2:
                # Patterns for heterocycles
                benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2n1")
                benzothiophene_pattern = Chem.MolFromSmarts("c1sc2ccccc2c1")

                # Check reactants for heterocycles
                heterocycle_count = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(benzimidazole_pattern) or mol.HasSubstructMatch(
                            benzothiophene_pattern
                        ):
                            heterocycle_count += 1

                # Check if product contains both heterocycles
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and heterocycle_count >= 2:
                    if product_mol.HasSubstructMatch(
                        benzimidazole_pattern
                    ) and product_mol.HasSubstructMatch(benzothiophene_pattern):
                        convergent_step_found = True
                        print("Found convergent heterocycle synthesis step")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_step_found
