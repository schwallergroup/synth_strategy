#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if the synthesis route includes formation of a tetracyclic structure.
    """

    def count_rings(mol):
        """Count the number of rings in a molecule."""
        if not mol:
            return 0
        ring_info = mol.GetRingInfo()
        return ring_info.NumRings()

    def has_tetracyclic_structure(node):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and count_rings(mol) >= 4:
                print(f"Tetracyclic structure detected: {node['smiles']}")
                return True

        # Check children
        for child in node.get("children", []):
            if has_tetracyclic_structure(child):
                return True

        return False

    def has_tetracyclic_formation(node):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                product_rings = count_rings(product_mol) if product_mol else 0

                # Check if any reactant has fewer rings than the product
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    reactant_rings = count_rings(reactant_mol) if reactant_mol else 0

                    if product_rings >= 4 and reactant_rings < product_rings:
                        print(f"Tetracyclic formation detected: {rsmi}")
                        return True
            except Exception as e:
                print(f"Error analyzing tetracyclic formation: {e}")

        # Check children
        for child in node.get("children", []):
            if has_tetracyclic_formation(child):
                return True

        return False

    # Check for either a tetracyclic structure or its formation
    return has_tetracyclic_structure(route) or has_tetracyclic_formation(route)
