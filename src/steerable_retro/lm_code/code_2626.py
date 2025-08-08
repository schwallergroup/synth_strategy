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
    Detects if the synthesis route employs a late-stage fragment coupling strategy
    where two complex fragments are joined in the final steps.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions (depth 0-1)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have at least 2 reactants (fragment coupling)
                if len(reactants) >= 2:
                    # Calculate complexity of reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    reactant_complexity = [mol.GetNumAtoms() for mol in reactant_mols if mol]

                    # Consider it a fragment coupling if both main reactants have significant complexity
                    if len(reactant_complexity) >= 2 and min(reactant_complexity[:2]) >= 8:
                        print(f"Found late-stage fragment coupling at depth {depth}")
                        result = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
