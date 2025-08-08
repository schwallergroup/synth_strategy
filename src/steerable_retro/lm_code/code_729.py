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
    This function detects a synthetic strategy involving both ring opening and
    ring formation steps in the same route.
    """
    has_ring_opening = False
    has_ring_formation = False

    def dfs_traverse(node):
        nonlocal has_ring_opening, has_ring_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Count rings in reactants and product
            reactant_rings = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    reactant_rings += len(rdmolops.GetSSSR(mol))

            product_mol = Chem.MolFromSmiles(product)
            product_rings = 0
            if product_mol:
                product_rings = len(rdmolops.GetSSSR(product_mol))

            print(f"Reaction: {rsmi}")
            print(f"Reactant rings: {reactant_rings}, Product rings: {product_rings}")

            # Check for ring opening or formation
            if product_rings < reactant_rings:
                has_ring_opening = True
                print(f"Detected ring opening: {rsmi}")
            elif product_rings > reactant_rings:
                has_ring_formation = True
                print(f"Detected ring formation: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Ring opening found: {has_ring_opening}, Ring formation found: {has_ring_formation}")
    return has_ring_opening and has_ring_formation
