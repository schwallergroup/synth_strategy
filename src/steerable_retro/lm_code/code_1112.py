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
    This function detects if the synthesis involves pyrazole ring formation from
    a phenylhydrazine and a diketone or similar precursors.
    """
    pyrazole_formed = False
    hydrazine_present = False
    diketone_present = False

    def dfs_traverse(node):
        nonlocal pyrazole_formed, hydrazine_present, diketone_present

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for pyrazole formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    pyrazole_patt = Chem.MolFromSmarts("c1nn[c]c1")
                    if product_mol.HasSubstructMatch(pyrazole_patt):
                        # Check if reactants contain hydrazine and diketone patterns
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                hydrazine_patt = Chem.MolFromSmarts("[NH2][NH]")
                                diketone_patt = Chem.MolFromSmarts(
                                    "[CX3](=[OX1])[CX4][CX3](=[OX1])"
                                )

                                if reactant_mol.HasSubstructMatch(hydrazine_patt):
                                    hydrazine_present = True
                                if reactant_mol.HasSubstructMatch(diketone_patt):
                                    diketone_present = True

                        if hydrazine_present and diketone_present:
                            pyrazole_formed = True
                            print("Detected pyrazole formation from hydrazine and diketone")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrazole_formed
