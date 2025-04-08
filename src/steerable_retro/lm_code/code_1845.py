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
    Detects if the synthesis includes a sequence of aromatic functionalizations
    (bromination, acylation, etc.) to establish substitution pattern.
    """
    aromatic_mods = []

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_mods

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic modifications
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Patterns for different aromatic modifications
                patterns = [
                    ("[c][Br]", "bromination"),
                    ("[c][Cl]", "chlorination"),
                    ("[c][C](=O)[C]", "acylation"),
                    ("[c][O][C]", "alkoxylation"),
                    ("[c][N]", "amination"),
                ]

                for pattern, mod_type in patterns:
                    patt = Chem.MolFromSmarts(pattern)
                    if product_mol.HasSubstructMatch(patt):
                        # Check if this is new (not in reactants)
                        mod_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(patt):
                                mod_in_reactants = True
                                break

                        if not mod_in_reactants:
                            aromatic_mods.append((depth, mod_type))
                            print(f"Aromatic {mod_type} detected at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort modifications by depth (early to late)
    aromatic_mods.sort(key=lambda x: x[0], reverse=True)

    # Check if we have at least 2 different aromatic modifications
    mod_types = set(mod[1] for mod in aromatic_mods)
    return len(mod_types) >= 2
