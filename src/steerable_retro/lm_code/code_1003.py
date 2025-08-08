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
    Detects if the synthetic route uses a linear fragment assembly strategy rather than convergent.
    """
    # Track the number of reactions with multiple complex reactants
    convergent_reactions = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal convergent_reactions, total_reactions

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            total_reactions += 1

            # Count complex reactants (those with more than 10 atoms)
            complex_reactants = 0
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 10:
                        complex_reactants += 1
                except:
                    print(f"Error processing reactant SMILES: {reactant}")

            # If there are multiple complex reactants, it's likely a convergent step
            if complex_reactants >= 2:
                convergent_reactions += 1
                print(f"Found convergent reaction step: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Calculate the ratio of convergent to total reactions
    if total_reactions > 0:
        convergent_ratio = convergent_reactions / total_reactions
        print(f"Convergent ratio: {convergent_ratio} ({convergent_reactions}/{total_reactions})")
        # If less than 30% of reactions are convergent, consider it a linear strategy
        return convergent_ratio < 0.3

    return False
