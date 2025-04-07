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
    Detects if the synthetic route contains a Friedel-Crafts acylation.
    """
    acylation_found = False

    def dfs_traverse(node):
        nonlocal acylation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Friedel-Crafts acylation
                acyl_chloride_pattern = Chem.MolFromSmarts("[CX3](=O)[Cl]")
                aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
                aryl_ketone_pattern = Chem.MolFromSmarts("c[C](=O)")

                has_acyl_chloride = False
                has_aromatic = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(acyl_chloride_pattern):
                            has_acyl_chloride = True
                        if reactant_mol.HasSubstructMatch(aromatic_pattern):
                            has_aromatic = True

                product_mol = Chem.MolFromSmiles(product)
                has_aryl_ketone = product_mol and product_mol.HasSubstructMatch(
                    aryl_ketone_pattern
                )

                if has_acyl_chloride and has_aromatic and has_aryl_ketone:
                    acylation_found = True
                    print("Found Friedel-Crafts acylation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return acylation_found
