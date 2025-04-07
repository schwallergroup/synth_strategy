#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the synthetic route involves an ester to hydrazide conversion.
    """
    ester_to_hydrazide_found = False

    def dfs_traverse(node):
        nonlocal ester_to_hydrazide_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester in reactants
                ester_pattern = Chem.MolFromSmarts("C(=O)OC")
                # Check for hydrazide in product
                hydrazide_pattern = Chem.MolFromSmarts("C(=O)NN")

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(ester_pattern):
                            prod_mol = Chem.MolFromSmiles(product)
                            if prod_mol and prod_mol.HasSubstructMatch(
                                hydrazide_pattern
                            ):
                                print("Ester to hydrazide conversion detected")
                                ester_to_hydrazide_found = True
                    except:
                        continue

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ester_to_hydrazide_found
