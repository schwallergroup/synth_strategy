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
    This function detects a synthetic strategy involving fragment coupling
    where two or more fragments are combined.
    """
    fragment_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Count number of reactant fragments
                reactants = reactants_part.split(".")
                num_reactants = len([r for r in reactants if r.strip()])

                # Count number of product fragments
                products = product_part.split(".")
                num_products = len([p for p in products if p.strip()])

                # If multiple reactants combine to fewer products, it's a coupling
                if num_reactants >= 2 and num_products < num_reactants:
                    print(
                        f"Found fragment coupling: {num_reactants} fragments → {num_products} fragments"
                    )
                    fragment_coupling_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return fragment_coupling_found
