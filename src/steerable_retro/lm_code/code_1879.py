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
    Detects if the synthesis route introduces a nitrile group via displacement of a leaving group.
    """
    has_nitrile_displacement = False

    def dfs_traverse(node):
        nonlocal has_nitrile_displacement

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for leaving group in reactants
            mesylate_pattern = Chem.MolFromSmarts("[#6]-OS(=O)(=O)C")
            tosylate_pattern = Chem.MolFromSmarts("[#6]-OS(=O)(=O)c1ccc(C)cc1")
            halide_pattern = Chem.MolFromSmarts("[#6]-[Cl,Br,I]")

            # Check for nitrile in product
            nitrile_pattern = Chem.MolFromSmarts("[#6]-C#N")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                has_leaving_group = any(
                    r.HasSubstructMatch(mesylate_pattern)
                    or r.HasSubstructMatch(tosylate_pattern)
                    or r.HasSubstructMatch(halide_pattern)
                    for r in reactant_mols
                    if r
                )
                has_nitrile = product_mol.HasSubstructMatch(nitrile_pattern)

                if has_leaving_group and has_nitrile:
                    print("Found nitrile introduction via leaving group displacement")
                    has_nitrile_displacement = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_nitrile_displacement
