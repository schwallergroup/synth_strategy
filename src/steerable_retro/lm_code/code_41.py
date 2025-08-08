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
    This function detects a linear synthesis strategy where each reaction
    has only one non-reagent reactant.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                # Count complex reactants (non-reagents)
                # Heuristic: complex molecules typically have aromatic rings or multiple functional groups
                complex_reactant_count = 0
                reactants = reactants_part.split(".")

                for reactant in reactants:
                    # Skip common reagents (simplified heuristic)
                    if len(reactant) > 10 and ("[c]" in reactant or "c1" in reactant):
                        complex_reactant_count += 1

                if complex_reactant_count > 1:
                    is_linear = False
                    print(
                        f"Found non-linear step with {complex_reactant_count} complex reactants: {rsmi}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return is_linear
