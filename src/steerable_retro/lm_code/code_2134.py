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
    This function detects if the synthetic route involves SNAr reactions
    forming C-O bonds (chloro-heterocycle + phenol/alcohol).
    """
    snar_reactions_count = 0

    def dfs_traverse(node):
        nonlocal snar_reactions_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr pattern: chloro-heterocycle + oxygen nucleophile
                chloro_hetero_pattern = Chem.MolFromSmarts("[n](:*):*[Cl]")
                oxygen_nucleophile_pattern = Chem.MolFromSmarts("[OX2H][#6]")

                chloro_hetero_present = False
                oxygen_nucleophile_present = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(chloro_hetero_pattern):
                        chloro_hetero_present = True

                    if reactant_mol.HasSubstructMatch(oxygen_nucleophile_pattern):
                        oxygen_nucleophile_present = True

                # Check if product has C-O bond where chlorine was
                if chloro_hetero_present and oxygen_nucleophile_present:
                    co_bond_pattern = Chem.MolFromSmarts("[n](:*):*[OX2][#6]")
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and product_mol.HasSubstructMatch(co_bond_pattern):
                        snar_reactions_count += 1
                        print(f"SNAr C-O bond formation detected (count: {snar_reactions_count})")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_reactions_count >= 1
