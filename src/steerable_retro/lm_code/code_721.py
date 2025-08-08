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
    Detects if the synthesis route is primarily linear but includes at least one
    convergent step where multiple complex fragments are combined.
    """
    has_convergent_step = False
    reaction_nodes = 0
    convergent_steps = 0

    def dfs_traverse(node):
        nonlocal has_convergent_step, reaction_nodes, convergent_steps

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            reaction_nodes += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if we have multiple complex reactants
            complex_reactants = 0

            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol:
                    # Consider both atom count and structural complexity
                    atom_count = mol.GetNumAtoms()
                    ring_count = 0
                    if mol.GetNumAtoms() > 3:  # Only calculate rings for non-trivial molecules
                        ring_info = mol.GetRingInfo()
                        if ring_info:
                            ring_count = ring_info.NumRings()

                    # Define complexity based on atoms and rings
                    if (atom_count > 8) or (atom_count > 5 and ring_count > 0):
                        complex_reactants += 1
                        print(
                            f"Complex reactant found: {r} with {atom_count} atoms and {ring_count} rings"
                        )

            if complex_reactants >= 2:
                # Try to extract depth, use a default if pattern not found
                depth = 999
                if "ID" in node["metadata"]:
                    depth_match = re.search(r"Depth: (\d+)", node["metadata"]["ID"])
                    if depth_match:
                        depth = int(depth_match.group(1))

                print(
                    f"Convergent step detected at depth {depth} with {complex_reactants} complex reactants"
                )
                has_convergent_step = True
                convergent_steps += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Route is primarily linear if it has at least one convergent step
    # but convergent steps are less than half of all reaction steps
    is_primarily_linear = has_convergent_step and (convergent_steps < reaction_nodes / 2)

    print(f"Total reaction nodes: {reaction_nodes}, Convergent steps: {convergent_steps}")
    print(f"Has convergent step: {has_convergent_step}, Is primarily linear: {is_primarily_linear}")

    return is_primarily_linear
