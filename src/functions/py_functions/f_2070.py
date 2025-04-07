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
    This function detects if the synthesis involves nucleophilic aromatic substitution
    with a chloropyridine and a nitrogen nucleophile.
    """
    nas_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nas_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check for chloropyridine pattern in reactants
                chloropyridine_pattern = Chem.MolFromSmarts(
                    "[#6]1:[#6]:[#6]:[#6](Cl):[#7]:[#6]:1"
                )

                # Check for amine pattern in reactants
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                # Check for N-substituted pyridine in product
                n_subst_pyridine_pattern = Chem.MolFromSmarts(
                    "[#6]1:[#6]:[#6]:[#6]([#7]):[#7]:[#6]:1"
                )

                if product_mol and product_mol.HasSubstructMatch(
                    n_subst_pyridine_pattern
                ):
                    chloro_found = False
                    amine_found = False

                    for r_mol in reactant_mols:
                        if r_mol:
                            if r_mol.HasSubstructMatch(chloropyridine_pattern):
                                chloro_found = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                amine_found = True

                    if chloro_found and amine_found:
                        print(
                            f"Found nucleophilic aromatic substitution at depth {depth}"
                        )
                        nas_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nas_detected
