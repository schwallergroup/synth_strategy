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
    This function detects a strategy where a halogen (Br, Cl, I) is installed early
    and then used in a late-stage coupling reaction.
    """
    early_halogenation = False
    halogen_used_in_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal early_halogenation, halogen_used_in_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for halogenation in early steps
                if depth >= 2:  # Early step
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            halogen_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,I]")
                            if product_mol.HasSubstructMatch(halogen_pattern):
                                print(f"Early halogenation detected at depth {depth}")
                                early_halogenation = True
                    except:
                        pass

                # Check for halogen used in coupling in late steps
                if depth <= 1:  # Late step
                    aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,I]")

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                                # Check if this is likely a coupling reaction
                                if any("B(O" in r for r in reactants) or any(
                                    "Sn(" in r for r in reactants
                                ):
                                    print(f"Halogen used in coupling at depth {depth}")
                                    halogen_used_in_coupling = True
                        except:
                            continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if halogen is installed early and used in late coupling
    return early_halogenation and halogen_used_in_coupling
