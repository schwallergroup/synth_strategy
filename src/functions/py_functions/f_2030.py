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
    This function detects Friedel-Crafts acylation on aromatic rings.
    It looks for the introduction of an acyl group onto an aromatic ring.
    """
    friedel_crafts_found = False

    def dfs_traverse(node):
        nonlocal friedel_crafts_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                try:
                    # Check if one reactant is an acyl chloride or similar
                    reactants = reactants_str.split(".")
                    acyl_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-Cl")
                    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")

                    has_acyl = False
                    has_aromatic = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(acyl_pattern):
                                has_acyl = True
                            if mol.HasSubstructMatch(aromatic_pattern):
                                has_aromatic = True

                    # Check if product has aromatic ring with acyl group
                    if has_acyl and has_aromatic:
                        product_mol = Chem.MolFromSmiles(product_str)
                        aromatic_acyl_pattern = Chem.MolFromSmarts("c-[#6](=[#8])-[#6]")

                        if product_mol and product_mol.HasSubstructMatch(
                            aromatic_acyl_pattern
                        ):
                            print("Friedel-Crafts acylation detected")
                            friedel_crafts_found = True
                except:
                    pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return friedel_crafts_found
