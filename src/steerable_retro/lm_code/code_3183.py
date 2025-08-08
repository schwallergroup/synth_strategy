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
    This function detects a synthetic strategy involving extension of a carbon chain
    (adding a significant carbon chain to a core structure).
    """
    chain_extension_found = False

    def dfs_traverse(node):
        nonlocal chain_extension_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a chain extension reaction
                # Look for patterns where a long alkyl chain (4+ carbons) is attached
                long_chain_patt = Chem.MolFromSmarts("[#6]-[#6]-[#6]-[#6]")

                product_mol = Chem.MolFromSmiles(product.split(".")[0])  # Take first product

                if product_mol and product_mol.HasSubstructMatch(long_chain_patt):
                    # Check if the chain is being attached to a core structure
                    # This is a simplification - in reality, you'd need to compare
                    # the product and reactants more carefully
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # If a reactant has an attachment point (like a halide)
                            attachment_patt = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")
                            if reactant_mol.HasSubstructMatch(attachment_patt):
                                print("Found chain extension reaction")
                                chain_extension_found = True
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return chain_extension_found
