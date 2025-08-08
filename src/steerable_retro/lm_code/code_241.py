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
    Detects if the synthesis follows a linear strategy rather than convergent
    (most reactions involve only one complex fragment)
    """
    # Track the number of reactions with multiple complex fragments
    convergent_reactions = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal convergent_reactions, total_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            total_reactions += 1

            # Count complex fragments (those with more than 10 atoms)
            complex_fragments = 0
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 10:
                        complex_fragments += 1
                except:
                    continue

            if complex_fragments >= 2:
                convergent_reactions += 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # If less than 25% of reactions are convergent, consider it a linear synthesis
    is_linear = total_reactions > 0 and (convergent_reactions / total_reactions) < 0.25
    if is_linear:
        print(
            f"Linear synthesis strategy detected ({convergent_reactions}/{total_reactions} convergent reactions)"
        )

    return is_linear
