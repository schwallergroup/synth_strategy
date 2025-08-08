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
    This function detects a convergent synthesis strategy where multiple complex fragments
    are combined throughout the route.
    """
    convergent_steps = 0

    def dfs_traverse(node):
        nonlocal convergent_steps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")
                print(f"Number of reactants: {len(reactants)}")

                # Check if multiple fragments are being joined
                if len(reactants) >= 2:
                    # Check complexity of fragments (having rings)
                    complex_fragments = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            num_rings = mol.GetRingInfo().NumRings()
                            print(f"Reactant {reactant} has {num_rings} rings")
                            if num_rings > 0:
                                complex_fragments += 1
                        else:
                            print(f"Failed to parse reactant SMILES: {reactant}")

                    print(f"Complex fragments: {complex_fragments}")
                    if complex_fragments >= 2:
                        print("Convergent step detected")
                        convergent_steps += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total convergent steps: {convergent_steps}")
    return convergent_steps >= 1
