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
    Detects if the synthesis involves alkylation of a phenol group,
    particularly with an epoxide.
    """
    phenol_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal phenol_alkylation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol pattern in reactants
                phenol_pattern = Chem.MolFromSmarts("c-[#8H]")

                # Check for epoxide pattern in reactants
                epoxide_pattern = Chem.MolFromSmarts("[#6]1[#8][#6]1")

                has_phenol = False
                has_epoxide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True
                            print(f"Found phenol in reactant at depth {depth}")
                        if mol.HasSubstructMatch(epoxide_pattern):
                            has_epoxide = True
                            print(f"Found epoxide in reactant at depth {depth}")

                # Check if product has aryl ether pattern
                product_mol = Chem.MolFromSmiles(product)
                aryl_ether_pattern = Chem.MolFromSmarts("c-[#8]-[#6]")

                if (
                    has_phenol
                    and has_epoxide
                    and product_mol
                    and product_mol.HasSubstructMatch(aryl_ether_pattern)
                ):
                    phenol_alkylation_detected = True
                    print(f"Phenol alkylation with epoxide detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return phenol_alkylation_detected
