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
    This function detects a strategy involving nucleophilic aromatic substitution
    where an aromatic chloride is displaced by an amine.
    """
    snar_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for aromatic chloride in reactants
                ar_cl_pattern = Chem.MolFromSmarts("[c]-[Cl]")

                # Check for C-N bond formation where C was previously connected to Cl
                # This is a simplified approach - a more robust implementation would track atom mappings
                reactants_mol = Chem.MolFromSmiles(reactants_part)
                product_mol = Chem.MolFromSmiles(product_part)

                if reactants_mol and product_mol:
                    if reactants_mol.HasSubstructMatch(ar_cl_pattern):
                        # Check if product has new C-N bond where C is aromatic
                        ar_n_pattern = Chem.MolFromSmarts("[c]-[#7]")
                        if product_mol.HasSubstructMatch(ar_n_pattern):
                            print(
                                f"Found potential nucleophilic aromatic substitution at depth {depth}"
                            )
                            snar_detected = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return snar_detected
