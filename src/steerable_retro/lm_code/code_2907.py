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
    Detects if the synthetic route involves amide bond formation between
    an acid/ester and an amine.
    """
    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester/acid and amine in reactants
                ester_or_acid_present = False
                amine_present = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Check for ester or acid
                        ester_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
                        if reactant_mol.HasSubstructMatch(ester_pattern):
                            ester_or_acid_present = True

                        # Check for amine
                        amine_pattern = Chem.MolFromSmarts("[NH2]")
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            amine_present = True

                # Check if product has amide bond
                if ester_or_acid_present and amine_present:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                        if product_mol.HasSubstructMatch(amide_pattern):
                            print(f"Found amide bond formation at depth {depth}")
                            found_amide_formation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_amide_formation
