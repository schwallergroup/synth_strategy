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
    Detects ketone protection as ketal/acetal in the synthesis route.
    """
    has_ketal_protection = False

    def dfs_traverse(node):
        nonlocal has_ketal_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ketal formation
                ketal_patt = Chem.MolFromSmarts("[#6]1[#8][#6][#6][#8]1")
                carbonyl_patt = Chem.MolFromSmarts("[#6]=[#8]")

                prod_mol = Chem.MolFromSmiles(product)

                # Check reactants for carbonyl
                has_carbonyl = False
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.HasSubstructMatch(carbonyl_patt):
                        has_carbonyl = True
                        break

                # Check if product has ketal
                if prod_mol and prod_mol.HasSubstructMatch(ketal_patt) and has_carbonyl:
                    print("Found ketone protection as ketal")
                    has_ketal_protection = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_ketal_protection
