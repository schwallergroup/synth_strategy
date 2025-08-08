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
    This function detects if the synthesis involves formation of a polycyclic system from a monocyclic aromatic starting material.
    """
    ring_formation_detected = False

    def dfs_traverse(node):
        nonlocal ring_formation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Count rings in reactants and product
            product_mol = Chem.MolFromSmiles(product)
            product_ring_count = 0
            if product_mol:
                product_ring_count = product_mol.GetRingInfo().NumRings()

            max_reactant_ring_count = 0
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    ring_count = reactant_mol.GetRingInfo().NumRings()
                    max_reactant_ring_count = max(max_reactant_ring_count, ring_count)

            # Check if product has more rings than any reactant
            if product_ring_count > max_reactant_ring_count:
                ring_formation_detected = True
                print(f"Ring formation detected: {rsmi}")
                print(
                    f"Product ring count: {product_ring_count}, Max reactant ring count: {max_reactant_ring_count}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Ring formation from aromatic strategy detected: {ring_formation_detected}")
    return ring_formation_detected
