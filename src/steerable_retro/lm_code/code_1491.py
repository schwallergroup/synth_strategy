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
    Detects if the synthesis involves halogenation (I, Br) to prepare a coupling partner.
    """
    found_halogenation = False

    def dfs_traverse(node):
        nonlocal found_halogenation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # Check if product contains aryl/heteroaryl halide
                product_mol = Chem.MolFromSmiles(product)
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I]")

                if product_mol and product_mol.HasSubstructMatch(aryl_halide_pattern):
                    # Check if reactants don't have the halide
                    has_halide_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                            has_halide_in_reactants = True
                            break

                    if not has_halide_in_reactants:
                        print("Found halogenation to prepare coupling partner")
                        found_halogenation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_halogenation
