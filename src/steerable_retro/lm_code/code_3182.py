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
    This function detects a synthetic strategy involving functionalization of aromatic systems.
    """
    aromatic_functionalization_found = False

    def dfs_traverse(node):
        nonlocal aromatic_functionalization_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this reaction involves modification of an aromatic system
                product_mol = Chem.MolFromSmiles(product.split(".")[0])  # Take first product

                if product_mol:
                    # Look for aromatic rings
                    aromatic_ring_patt = Chem.MolFromSmarts("c1ccccc1")

                    if product_mol.HasSubstructMatch(aromatic_ring_patt):
                        # Check if any functional group is attached to the aromatic ring
                        # Common functional groups attached to aromatic rings
                        functional_groups = [
                            "[c]-[#7]",  # Aromatic amine
                            "[c]-[#16]",  # Aromatic sulfur
                            "[c]-[#8]",  # Aromatic oxygen
                            "[c]-[#6](=[#8])",  # Aromatic carbonyl
                        ]

                        for fg in functional_groups:
                            fg_patt = Chem.MolFromSmarts(fg)
                            if product_mol.HasSubstructMatch(fg_patt):
                                print(f"Found aromatic functionalization with pattern {fg}")
                                aromatic_functionalization_found = True
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return aromatic_functionalization_found
