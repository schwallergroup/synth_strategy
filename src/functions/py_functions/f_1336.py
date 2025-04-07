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
    Detects a convergent synthesis strategy using Friedel-Crafts acylation
    to combine two complex fragments, followed by ketone reduction.
    """
    # Track if we found the key reactions
    found_friedel_crafts = False
    found_ketone_reduction = False

    def dfs_traverse(node):
        nonlocal found_friedel_crafts, found_ketone_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Friedel-Crafts acylation
                if len(reactants) == 2:
                    acid_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")
                    aromatic_pattern = Chem.MolFromSmarts("[c]")
                    ketone_pattern = Chem.MolFromSmarts("[c][C](=O)[c]")

                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        product_mol
                        and any(
                            m and m.HasSubstructMatch(acid_chloride_pattern)
                            for m in reactant_mols
                        )
                        and any(
                            m and m.HasSubstructMatch(aromatic_pattern)
                            for m in reactant_mols
                        )
                        and product_mol.HasSubstructMatch(ketone_pattern)
                    ):
                        print("Found Friedel-Crafts acylation")
                        found_friedel_crafts = True

                # Check for ketone reduction
                ketone_pattern = Chem.MolFromSmarts("[c][C](=O)[c]")
                methylene_pattern = Chem.MolFromSmarts("[c][CH2][c]")

                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # In forward direction, this would be a ketone reduction
                # In retrosynthetic direction (as shown), this is an oxidation
                if (
                    product_mol
                    and any(
                        m and m.HasSubstructMatch(ketone_pattern) for m in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(methylene_pattern)
                ):
                    print("Found ketone reduction (in forward direction)")
                    found_ketone_reduction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both key reactions are found
    return found_friedel_crafts and found_ketone_reduction
