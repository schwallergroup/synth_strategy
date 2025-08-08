#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if the synthesis route uses a late-stage Negishi coupling to introduce a cyclobutyl group.
    Late stage means within the first 2 steps (depth 0-1).
    """
    found_negishi = False

    def dfs_traverse(node, depth=0):
        nonlocal found_negishi

        if node["type"] == "reaction" and depth <= 1:
            # Check if this is a Negishi coupling reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for zinc reagent
                zinc_present = any("[Zn+]" in reactant for reactant in reactants)

                # Check for cyclobutyl group
                cyclobutyl_pattern = Chem.MolFromSmarts("[CH]1[CH2][CH2][CH2]1")
                cyclobutyl_present = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(cyclobutyl_pattern):
                            cyclobutyl_present = True
                            break
                    except:
                        continue

                if zinc_present and cyclobutyl_present:
                    print(f"Found Negishi coupling with cyclobutyl at depth {depth}")
                    found_negishi = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_negishi
