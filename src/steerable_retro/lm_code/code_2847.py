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
    This function detects Suzuki coupling for biaryl formation.
    """
    boronic_acid_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
    aryl_halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#35,#53,#17]")
    biaryl_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#6]:[#6]")
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for boronic acid and aryl halide in reactants
                has_boronic_acid = False
                has_aryl_halide = False

                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True

                    if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

                # Check for biaryl in product
                product_mol = Chem.MolFromSmiles(product_smiles)
                has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)

                if has_boronic_acid and has_aryl_halide and has_biaryl:
                    print("Suzuki coupling detected")
                    suzuki_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
