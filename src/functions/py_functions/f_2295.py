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
    This function detects ether formation between aromatic ring and cyclic ether.
    """
    ether_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains aromatic ether with cyclic component
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("c-[O]-C1CCOC1")
                ):
                    # Check if reactants contain phenol and cyclic component
                    has_phenol = False
                    has_cyclic_component = False

                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("c[OH]")
                            ):
                                has_phenol = True
                            if reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("C1CCOC1")
                            ):
                                has_cyclic_component = True

                    if has_phenol and has_cyclic_component:
                        ether_formation_detected = True
                        print(f"Ether formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ether_formation_detected
