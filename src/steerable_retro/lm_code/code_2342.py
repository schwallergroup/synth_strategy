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
    This function detects if the route contains a fragment coupling strategy
    where two or more complex fragments are combined.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Count number of fragments in reactants
            reactant_fragments = reactants_part.split(".")

            # Only consider reactions with multiple reactant fragments
            if len(reactant_fragments) >= 2:
                # Check complexity of fragments (simplified: count atoms)
                complex_fragments = 0
                for frag in reactant_fragments:
                    # Count non-hydrogen atoms as a simple complexity measure
                    atom_count = sum(1 for c in frag if c.isupper())
                    if atom_count >= 5:  # Consider fragments with 5+ atoms as complex
                        complex_fragments += 1

                if complex_fragments >= 2:
                    print(f"Found fragment coupling with {complex_fragments} complex fragments")
                    found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return found
