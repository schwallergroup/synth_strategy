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
    This function detects a linear synthesis with late-stage heterocycle incorporation.
    Checks if the synthesis is linear (no convergent steps) and if a heterocycle
    (like piperazine) is introduced in the final steps.
    """
    is_linear = True
    late_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, late_heterocycle

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if more than 2 significant reactants (would indicate convergent synthesis)
                significant_reactants = 0
                for reactant in reactants:
                    if len(reactant) > 5:  # Simple heuristic to filter out small reagents
                        significant_reactants += 1

                if significant_reactants > 2:
                    is_linear = False

                # Check for heterocycle in late-stage reactions
                if depth <= 1:  # Only consider late-stage reactions (depth 0 or 1)
                    piperazine_pattern = Chem.MolFromSmarts("N1CCN(C)CC1")

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(piperazine_pattern):
                                late_heterocycle = True
                                print("Detected late-stage heterocycle incorporation")
                                break
                        except:
                            continue

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear and late_heterocycle
