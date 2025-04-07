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
    This function detects if the synthesis route includes conversion of an alkyl chloride to an alkene.
    """
    # Track if we found the conversion
    found_conversion = False

    def dfs_traverse(node, depth=0):
        nonlocal found_conversion

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alkyl chloride to alkene conversion
                alkyl_chloride_pattern = Chem.MolFromSmarts("[#6]-[#17]")
                alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")

                # Check reactant for alkyl chloride
                has_alkyl_chloride = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(alkyl_chloride_pattern):
                            has_alkyl_chloride = True
                    except:
                        continue

                # Check product for alkene
                has_alkene = False
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(alkene_pattern):
                        has_alkene = True
                except:
                    pass

                if has_alkyl_chloride and has_alkene:
                    found_conversion = True
                    print("Found alkyl chloride to alkene conversion")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_conversion
