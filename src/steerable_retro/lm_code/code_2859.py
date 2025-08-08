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
    Detects if the synthesis contains multiple Suzuki coupling reactions
    """
    suzuki_coupling_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_count

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for Suzuki coupling patterns
                # Look for aryl bromide and boronic acid/ester derivatives
                if "Br" in reactants_str and (
                    "B(O" in reactants_str or "BO" in reactants_str or "OB" in reactants_str
                ):

                    reactants = reactants_str.split(".")

                    # Try to identify aryl bromide and boronic acid/ester
                    aryl_bromide = None
                    boronic_derivative = None

                    for reactant in reactants:
                        if "Br" in reactant:
                            aryl_bromide = Chem.MolFromSmiles(reactant)
                        if "B(O" in reactant or "BO" in reactant or "OB" in reactant:
                            boronic_derivative = Chem.MolFromSmiles(reactant)

                    product_mol = Chem.MolFromSmiles(product)

                    if aryl_bromide and boronic_derivative and product_mol:
                        # Check if product has more C-C bonds than reactants
                        suzuki_coupling_count += 1
                        print(f"Found Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found at least 2 Suzuki couplings
    return suzuki_coupling_count >= 2
