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
    Detects if the synthesis uses a convergent approach where two complex fragments
    are prepared separately and then connected in a late-stage coupling reaction.
    """
    # Initialize tracking variables
    fragment_count = 0
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal fragment_count, late_stage_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for late-stage coupling (depth 0 or 1)
            if depth <= 1:
                # Count complex reactants (those with more than 10 atoms)
                complex_reactants = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > 10:
                            complex_reactants += 1
                    except:
                        continue

                # If we have at least 2 complex reactants, it's likely a convergent coupling
                if complex_reactants >= 2:
                    fragment_count = complex_reactants
                    late_stage_coupling = True
                    print(
                        f"Detected late-stage coupling with {complex_reactants} complex fragments at depth {depth}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if we have a late-stage coupling with at least 2 fragments
    strategy_present = late_stage_coupling and fragment_count >= 2
    print(f"Convergent synthesis strategy detected: {strategy_present}")
    return strategy_present
