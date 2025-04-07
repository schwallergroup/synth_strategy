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
    Detects a strategy involving N-alkylation of a lactam, specifically
    the reaction of a bromomethyl compound with a lactam to form an N-alkylated product.
    """
    n_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")

            # Extract reactants and product
            parts = rsmi.split(">")
            if len(parts) >= 3:
                reactants = parts[0].split(".")
                product = parts[2]

                # Check for bromomethyl compound
                bromomethyl_pattern = Chem.MolFromSmarts("[c]-[CH2][Br]")

                # Check for lactam pattern
                lactam_pattern = Chem.MolFromSmarts("[NH]-[C](=[O])-[C]")

                # Check for N-alkylated lactam pattern
                n_alkylated_pattern = Chem.MolFromSmarts("[N](-[CH2][c])-[C](=[O])-[C]")

                has_bromomethyl = False
                has_lactam = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(bromomethyl_pattern):
                            has_bromomethyl = True
                        if reactant_mol.HasSubstructMatch(lactam_pattern):
                            has_lactam = True

                product_mol = Chem.MolFromSmiles(product)
                has_n_alkylated = False
                if product_mol:
                    has_n_alkylated = product_mol.HasSubstructMatch(n_alkylated_pattern)

                # If all conditions are met, it's an N-alkylation of lactam
                if has_bromomethyl and has_lactam and has_n_alkylated:
                    print("Detected N-alkylation of lactam")
                    n_alkylation_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return n_alkylation_detected
