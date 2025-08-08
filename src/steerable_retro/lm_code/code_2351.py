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
    Detects a synthetic strategy involving halogenation of aromatic rings.
    """
    halogenation_detected = False

    def dfs_traverse(node):
        nonlocal halogenation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for halogenation pattern
            aryl_halide_pattern = Chem.MolFromSmarts("c1ccccc1[I,Br,Cl]")

            # Check if product contains aryl halide but reactants don't
            if (
                product
                and Chem.MolFromSmiles(product)
                and Chem.MolFromSmiles(product).HasSubstructMatch(aryl_halide_pattern)
            ):
                reactants_have_halide = False
                for reactant in reactants:
                    if Chem.MolFromSmiles(reactant) and Chem.MolFromSmiles(
                        reactant
                    ).HasSubstructMatch(aryl_halide_pattern):
                        reactants_have_halide = True
                        break

                if not reactants_have_halide:
                    halogenation_detected = True
                    print("Aromatic halogenation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return halogenation_detected
