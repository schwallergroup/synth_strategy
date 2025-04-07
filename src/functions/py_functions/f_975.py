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
    This function detects if the synthesis uses an aryl ether linker installation
    strategy via phenol alkylation.
    """
    has_aryl_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aryl_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol alkylation
                phenol_pattern = Chem.MolFromSmarts("c[O;H1]")
                bromo_pattern = Chem.MolFromSmarts("[Br][C]")
                aryl_ether_pattern = Chem.MolFromSmarts("c[O][C][C][C][C]C(=O)")

                has_phenol = False
                has_bromo = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True
                        if mol.HasSubstructMatch(bromo_pattern):
                            has_bromo = True

                product_mol = Chem.MolFromSmiles(product)
                has_aryl_ether = product_mol and product_mol.HasSubstructMatch(
                    aryl_ether_pattern
                )

                if has_phenol and has_bromo and has_aryl_ether:
                    print(
                        f"Detected aryl ether formation via phenol alkylation at depth {depth}"
                    )
                    has_aryl_ether_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_aryl_ether_formation
