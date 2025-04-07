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
    Detects if the synthesis route involves late-stage introduction of a tertiary amine-containing group
    (like dimethylaminoethyl) via alkylation of a hydroxyl group in the final 1-2 steps.
    """
    tertiary_amine_introduced = False
    hydroxyl_alkylation = False
    late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal tertiary_amine_introduced, hydroxyl_alkylation, late_stage

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions in the final 2 steps
            late_stage = True

            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for tertiary amine pattern in reactants
                tertiary_amine_pattern = Chem.MolFromSmarts("[#6][#6][#7]([#6])[#6]")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(tertiary_amine_pattern):
                            tertiary_amine_introduced = True
                            print(f"Found tertiary amine in reactant: {reactant}")
                            break
                    except:
                        continue

                # Check for hydroxyl alkylation
                try:
                    for child in node.get("children", []):
                        if child["type"] == "mol":
                            reactant_mol = Chem.MolFromSmiles(child["smiles"])
                            if reactant_mol:
                                hydroxyl_pattern = Chem.MolFromSmarts("[#8H]")
                                if reactant_mol.HasSubstructMatch(hydroxyl_pattern):
                                    product_mol = Chem.MolFromSmiles(product)
                                    ether_pattern = Chem.MolFromSmarts("[#8][#6][#6]")
                                    if product_mol and product_mol.HasSubstructMatch(
                                        ether_pattern
                                    ):
                                        hydroxyl_alkylation = True
                                        print(
                                            f"Found hydroxyl alkylation: {child['smiles']} -> {product}"
                                        )
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = tertiary_amine_introduced and hydroxyl_alkylation and late_stage
    print(f"Late-stage tertiary amine introduction detected: {result}")
    return result
