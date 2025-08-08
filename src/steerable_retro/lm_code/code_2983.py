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
    Detects if the synthesis route involves a malonate alkylation strategy
    (use of diethyl malonate as a building block for C-C bond formation).
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if any reactant is diethyl malonate or similar
                for reactant in reactants:
                    try:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        # Pattern for diethyl malonate or similar malonate derivatives
                        malonate_pattern = Chem.MolFromSmarts("C(=O)OCC.C(=O)OCC")
                        if reactant_mol and reactant_mol.HasSubstructMatch(malonate_pattern):
                            print(f"Detected malonate derivative at depth {depth}")
                            result = True
                            break
                    except:
                        continue

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
