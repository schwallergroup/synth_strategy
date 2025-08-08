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
    This function detects a linear synthesis strategy building a complex
    multi-aromatic system with ether linkages.
    """
    aromatic_ring_count = 0
    ether_linkage_count = 0
    linear_structure = True

    def count_aromatic_rings(mol):
        """Count the number of aromatic rings in a molecule"""
        if not mol:
            return 0

        # Find all aromatic rings
        ri = mol.GetRingInfo()
        aromatic_rings = 0

        for ring in ri.AtomRings():
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_rings += 1

        return aromatic_rings

    def dfs_traverse(node):
        nonlocal aromatic_ring_count, ether_linkage_count, linear_structure

        if node["type"] == "mol" and node.get("in_stock", False):
            # This is a starting material
            mol = Chem.MolFromSmiles(node["smiles"])
            aromatic_ring_count += count_aromatic_rings(mol)

        if node["type"] == "reaction":
            # Check if this is a branching reaction (more than one reactant that's not in stock)
            non_stock_reactants = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_stock_reactants += 1

            if non_stock_reactants > 1:
                linear_structure = False

            # Check for ether linkages in the product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product_mol = Chem.MolFromSmiles(product_smiles)

            ether_pattern = Chem.MolFromSmarts("cOc")
            if product_mol and product_mol.HasSubstructMatch(ether_pattern):
                matches = product_mol.GetSubstructMatches(ether_pattern)
                ether_linkage_count += len(matches)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if it's a linear synthesis with multiple aromatic rings and ether linkages
    return linear_structure and aromatic_ring_count >= 2 and ether_linkage_count >= 1
