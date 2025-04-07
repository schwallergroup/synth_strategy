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
    This function detects if the synthetic route involves a specific sequence of
    aromatic functional group transformations: nitro → amine → halogen → biaryl.
    """
    # Track if we've seen each transformation in the correct order
    transformations_seen = []

    def dfs_traverse(node):
        nonlocal transformations_seen

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            try:
                reactant_mol = Chem.MolFromSmiles(reactants_part)
                product_mol = Chem.MolFromSmiles(products_part)

                if reactant_mol and product_mol:
                    # Check for nitro to amine transformation
                    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
                    halogen_pattern = Chem.MolFromSmarts("[#6]-[#9,#17,#35,#53]")

                    # Detect nitro to amine
                    if reactant_mol.HasSubstructMatch(
                        nitro_pattern
                    ) and product_mol.HasSubstructMatch(amine_pattern):
                        transformations_seen.append("nitro_to_amine")
                        print("Nitro to amine transformation detected")

                    # Detect amine to halogen
                    elif reactant_mol.HasSubstructMatch(
                        amine_pattern
                    ) and product_mol.HasSubstructMatch(halogen_pattern):
                        transformations_seen.append("amine_to_halogen")
                        print("Amine to halogen transformation detected")

                    # Detect biaryl coupling (simplified check)
                    elif len(reactants_part.split(".")) >= 2:
                        # If we have multiple reactants and one product, it might be a coupling
                        if "." not in products_part:
                            transformations_seen.append("biaryl_coupling")
                            print("Potential biaryl coupling detected")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if we have the correct sequence
    correct_sequence = ["biaryl_coupling", "amine_to_halogen", "nitro_to_amine"]

    # Check if all elements of correct_sequence appear in transformations_seen in the right order
    last_index = -1
    for item in correct_sequence:
        if item not in transformations_seen:
            return False
        current_index = transformations_seen.index(item)
        if current_index <= last_index:
            return False
        last_index = current_index

    return True
