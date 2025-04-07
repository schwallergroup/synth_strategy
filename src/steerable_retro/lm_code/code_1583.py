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
    Detects if the synthetic route involves formation of an aryl ether.
    """
    has_aryl_ether_formation = False

    def dfs_traverse(node):
        nonlocal has_aryl_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol in reactants
                phenol_pattern = Chem.MolFromSmarts("[c]-[OH]")

                has_phenol = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(phenol_pattern):
                        has_phenol = True
                        break

                # Check for aryl ether in product
                if has_phenol:
                    aryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[C]")
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(aryl_ether_pattern):
                        print("Detected aryl ether formation")
                        has_aryl_ether_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_aryl_ether_formation
