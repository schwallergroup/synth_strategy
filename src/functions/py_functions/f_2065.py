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
    This function detects if the synthesis route involves C-C bond formation
    using organometallic reagents (Grignard or organotin).
    """
    has_organometallic_cc = False

    def dfs_traverse(node):
        nonlocal has_organometallic_cc

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for Grignard reagents
                grignard_pattern = Chem.MolFromSmarts("[#6]-[Mg]")

                # Check for organotin reagents
                organotin_pattern = Chem.MolFromSmarts("[#6]-[Sn]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(grignard_pattern):
                            has_organometallic_cc = True
                            print("Grignard reagent detected in synthesis")
                        elif reactant_mol.HasSubstructMatch(organotin_pattern):
                            has_organometallic_cc = True
                            print("Organotin reagent detected in synthesis")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_organometallic_cc
