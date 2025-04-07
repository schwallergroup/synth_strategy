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
    This function detects if the synthetic route involves aromatic halogenation
    in the presence of a nitro group.
    """
    # Track if we found aromatic halogenation and nitro group
    found_halogenation = False
    has_nitro_group = False

    def dfs_traverse(node, depth=0):
        nonlocal found_halogenation, has_nitro_group

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for halogenation (C-H to C-Cl/Br/F/I)
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    halogen_pattern = Chem.MolFromSmarts("c-[Cl,Br,F,I]")
                    if product_mol.HasSubstructMatch(halogen_pattern):
                        # Check if this is a new halogen (not in reactants)
                        halogen_count_product = len(
                            product_mol.GetSubstructMatches(halogen_pattern)
                        )

                        total_halogen_count_reactants = 0
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if (
                                reactant_mol
                                and "Cl" in reactant
                                or "Br" in reactant
                                or "F" in reactant
                                or "I" in reactant
                            ):
                                total_halogen_count_reactants += len(
                                    reactant_mol.GetSubstructMatches(halogen_pattern)
                                )

                        if halogen_count_product > total_halogen_count_reactants:
                            print(f"Found halogenation at depth {depth}")
                            found_halogenation = True

                # Check for nitro group in any molecule
                nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                    has_nitro_group = True
                    print(f"Found nitro group at depth {depth}")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                        has_nitro_group = True
                        print(f"Found nitro group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return found_halogenation and has_nitro_group
