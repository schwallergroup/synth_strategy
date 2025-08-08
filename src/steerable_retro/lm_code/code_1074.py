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
    Detects if the synthesis route follows a linear strategy rather than convergent.
    A linear strategy typically has only one complex reactant per step.
    """
    # Count the number of steps with multiple complex reactants
    convergent_steps = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal convergent_steps, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants (more than 10 atoms)
                complex_reactant_count = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 10:
                        complex_reactant_count += 1

                if complex_reactant_count >= 2:
                    convergent_steps += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If less than 20% of steps are convergent, consider it a linear synthesis
    is_linear = total_steps > 0 and (convergent_steps / total_steps) < 0.2
    if is_linear:
        print(
            f"Linear synthesis detected: {convergent_steps} convergent steps out of {total_steps} total steps"
        )

    return is_linear
