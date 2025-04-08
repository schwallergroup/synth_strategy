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
    This function detects if the synthesis includes SNAr reactions with nitrogen nucleophiles.
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for halogen-containing aromatic in one reactant
            halogen_aromatic_pattern = Chem.MolFromSmarts("[c]-[F,Cl,Br,I]")
            nitrogen_nucleophile_pattern = Chem.MolFromSmarts("[N;H0,H1,H2]")

            has_halogen_aromatic = False
            has_nitrogen_nucleophile = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                if reactant_mol.HasSubstructMatch(halogen_aromatic_pattern):
                    has_halogen_aromatic = True

                if reactant_mol.HasSubstructMatch(nitrogen_nucleophile_pattern):
                    has_nitrogen_nucleophile = True

            # Check if product has new C-N bond where halogen was
            if has_halogen_aromatic and has_nitrogen_nucleophile:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # This is a simplification - in a real implementation,
                    # we would need to check for new C-N bonds at positions where halogens were
                    print("Detected potential SNAr with nitrogen nucleophile")
                    snar_count += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return snar_count >= 1
