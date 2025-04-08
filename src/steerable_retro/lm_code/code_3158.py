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
    Detects a synthesis strategy involving nucleophilic aromatic substitution with a thiol.
    """
    has_snar_with_thiol = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar_with_thiol

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product = Chem.MolFromSmiles(product_smiles)

                # Check for SNAr with thiol
                chloro_aromatic_pattern = Chem.MolFromSmarts("[c:1][Cl:2]")
                thiol_pattern = Chem.MolFromSmarts("[SH]")
                cs_bond_pattern = Chem.MolFromSmarts("[c:1][S:2]")

                has_chloro_aromatic = False
                has_thiol = False

                for reactant in reactants:
                    if reactant:
                        if reactant.HasSubstructMatch(chloro_aromatic_pattern):
                            has_chloro_aromatic = True
                        if reactant.HasSubstructMatch(thiol_pattern):
                            has_thiol = True

                if (
                    has_chloro_aromatic
                    and has_thiol
                    and product
                    and product.HasSubstructMatch(cs_bond_pattern)
                ):
                    has_snar_with_thiol = True
                    print("Found SNAr with thiol")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"SNAr with thiol strategy: {has_snar_with_thiol}")
    return has_snar_with_thiol
