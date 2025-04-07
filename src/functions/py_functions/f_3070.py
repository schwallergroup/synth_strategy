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
    Detects N-arylation coupling between a halogenated heterocycle and an NH-heterocycle.
    """
    found_n_arylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_arylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-arylation pattern
                nh_heterocycle_pattern = Chem.MolFromSmarts("[nH]")
                halogen_heterocycle_pattern = Chem.MolFromSmarts("c-[Br,Cl,I,F]")
                n_arylated_pattern = Chem.MolFromSmarts("c-n")

                # Check reactants for NH-heterocycle and halogenated heterocycle
                has_nh_heterocycle = False
                has_halogen_heterocycle = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(nh_heterocycle_pattern):
                            has_nh_heterocycle = True
                        if reactant_mol.HasSubstructMatch(halogen_heterocycle_pattern):
                            has_halogen_heterocycle = True

                # Check product for N-arylated heterocycle
                product_mol = Chem.MolFromSmiles(product)
                has_n_arylated = product_mol and product_mol.HasSubstructMatch(
                    n_arylated_pattern
                )

                if has_nh_heterocycle and has_halogen_heterocycle and has_n_arylated:
                    found_n_arylation = True
                    print(f"Found N-arylation coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_n_arylation
