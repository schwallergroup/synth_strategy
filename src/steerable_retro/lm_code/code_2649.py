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
    This function detects if the synthetic route follows a linear synthesis strategy
    (no convergent steps with multiple fragments combining).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If there are multiple reactants (excluding reagents), it's potentially convergent
                if len(reactants) > 1:
                    # Check if these are actual fragments or just reagents
                    significant_reactants = 0
                    for reactant in reactants:
                        # Simple heuristic: reactants with more than 10 atoms are considered significant
                        # This could be refined based on your specific needs
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.GetNumAtoms() > 10:
                                significant_reactants += 1
                        except:
                            print(f"Error processing reactant SMILES: {reactant}")

                    if significant_reactants > 1:
                        print(f"Found convergent step with multiple significant reactants: {rsmi}")
                        is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
