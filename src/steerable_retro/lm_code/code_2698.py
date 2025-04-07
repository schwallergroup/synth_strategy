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
    Detects the use of Grignard reagents in the synthesis.
    """
    has_grignard_reaction = False

    def dfs_traverse(node):
        nonlocal has_grignard_reaction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for Grignard reagent pattern
                grignard_patt = Chem.MolFromSmarts("[#6][Mg][Br,Cl,I,F]")

                for r in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(grignard_patt):
                            print("Found Grignard reagent in reaction")
                            has_grignard_reaction = True
                            break
                    except:
                        # Some Grignard reagents might not parse well
                        if "[Mg]" in r:
                            print("Found potential Grignard reagent by text search")
                            has_grignard_reaction = True
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_grignard_reaction
