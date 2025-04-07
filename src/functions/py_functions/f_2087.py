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
    Detects if a nitrile group is present throughout the synthesis but not modified.
    """
    nitrile_present = False
    nitrile_modified = False

    def dfs_traverse(node):
        nonlocal nitrile_present, nitrile_modified

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Check for nitrile pattern
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                if mol.HasSubstructMatch(nitrile_pattern):
                    nitrile_present = True

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if nitrile is present in reactants but not in product
            nitrile_in_reactants = False
            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                        if mol.HasSubstructMatch(nitrile_pattern):
                            nitrile_in_reactants = True

            prod_mol = Chem.MolFromSmiles(product)
            nitrile_in_product = False
            if prod_mol:
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                if prod_mol.HasSubstructMatch(nitrile_pattern):
                    nitrile_in_product = True

            if nitrile_in_reactants != nitrile_in_product:
                nitrile_modified = True
                print("Nitrile group was modified in a reaction")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitrile_present and not nitrile_modified
