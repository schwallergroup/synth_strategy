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
    This function detects SNAr reaction of halide with amine in the synthesis route.
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aromatic halide in reactants
                aromatic_halide_pattern = Chem.MolFromSmarts("[c]-[F,Cl,Br,I]")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                # Check if reactants have aromatic halide and amine
                has_aromatic_halide = False
                has_amine = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if aromatic_halide_pattern and reactant_mol.HasSubstructMatch(
                        aromatic_halide_pattern
                    ):
                        has_aromatic_halide = True

                    if amine_pattern and reactant_mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

                # Check if product has new C-N bond where halide was
                if has_aromatic_halide and has_amine:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Look for aromatic carbon connected to nitrogen
                        c_n_pattern = Chem.MolFromSmarts("[c]-[#7;!$(N~[!#6])]")
                        if c_n_pattern and product_mol.HasSubstructMatch(c_n_pattern):
                            snar_found = True
                            print("SNAr reaction detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_found
