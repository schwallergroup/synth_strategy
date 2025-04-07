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
    This function detects if the synthetic route involves a late-stage SNAr reaction.
    Looks for aryl halide and amine coupling in the final steps of the synthesis.
    """
    snar_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_detected

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions at depth 0 or 1 (late stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for aryl halide in reactants
                reactants = Chem.MolFromSmiles(reactants_smiles)
                if reactants:
                    aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Br,Cl,F,I]")
                    amine_pattern = Chem.MolFromSmarts("[NH]")

                    if reactants.HasSubstructMatch(
                        aryl_halide_pattern
                    ) and reactants.HasSubstructMatch(amine_pattern):

                        # Check for C-N bond formation in product
                        product = Chem.MolFromSmiles(product_smiles)
                        if product:
                            c_n_bond_pattern = Chem.MolFromSmarts("[c]-[N]")
                            if product.HasSubstructMatch(c_n_bond_pattern):
                                snar_detected = True
                                print(f"Detected late-stage SNAr reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Late-stage SNAr strategy detected: {snar_detected}")
    return snar_detected
