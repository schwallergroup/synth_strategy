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
    This function detects if the synthesis involves building or modifying a hydantoin scaffold.
    """
    # Track if we found the hydantoin core
    found_hydantoin = False
    # Track if hydantoin is present in multiple steps
    hydantoin_step_count = 0

    def dfs_traverse(node):
        nonlocal found_hydantoin, hydantoin_step_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for hydantoin pattern in product
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        hydantoin_pattern = Chem.MolFromSmarts(
                            "[#7]-1-[#6](=[#8])-[#7]-[#6](=[#8])-[#6]-1"
                        )
                        if prod_mol.HasSubstructMatch(hydantoin_pattern):
                            found_hydantoin = True
                            hydantoin_step_count += 1
                            print(f"Found hydantoin in step: {rsmi}")
                except:
                    pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found hydantoin in multiple steps (indicating it's a key scaffold)
    result = found_hydantoin and hydantoin_step_count >= 2
    print(
        f"Hydantoin scaffold strategy detected: {result} (hydantoin_step_count={hydantoin_step_count})"
    )
    return result
