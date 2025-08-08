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
    This function detects if the synthesis involves sequential functionalization
    of an aromatic ring at different positions.
    """
    # Track functionalization events
    functionalization_events = []

    def dfs_traverse(node, depth=0):
        nonlocal functionalization_events

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check for various aromatic functionalizations
                    functional_groups = [
                        ("[c]-[N+](=[O])[O-]", "nitro"),
                        ("[c]-[NH2]", "amine"),
                        ("[c]-[NH]-[C](=[O])-[CH3]", "acetamide"),
                        ("[c]-[OH]", "hydroxyl"),
                        ("[c]-[O]-[CH3]", "methoxy"),
                        ("[c]-[C](=[O])-[CH3]", "acetyl"),
                        ("[c]-[Br]", "bromo"),
                        ("[c]-[Cl]", "chloro"),
                        ("[c]-[I]", "iodo"),
                    ]

                    for smarts, name in functional_groups:
                        # Check if a functional group is added
                        if not reactants_mol.HasSubstructMatch(
                            Chem.MolFromSmarts(smarts)
                        ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                            functionalization_events.append((depth, name, "addition"))
                            print(f"Detected {name} addition at depth {depth}")

                        # Check if a functional group is removed
                        if reactants_mol.HasSubstructMatch(
                            Chem.MolFromSmarts(smarts)
                        ) and not product_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                            functionalization_events.append((depth, name, "removal"))
                            print(f"Detected {name} removal at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort events by depth
    functionalization_events.sort(key=lambda x: x[0])

    # Check if we have at least 3 different functionalization events
    unique_events = set([(name, action) for _, name, action in functionalization_events])
    sequential_functionalization = len(unique_events) >= 3

    if sequential_functionalization:
        print(
            f"Sequential aromatic functionalization detected with {len(unique_events)} different events"
        )

    return sequential_functionalization
