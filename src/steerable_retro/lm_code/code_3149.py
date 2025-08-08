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
    Detects a linear synthesis strategy building a piperazine scaffold with a
    functionalized benzyl side chain.
    """
    found_piperazine = False
    found_fluorobenzyl = False
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal found_piperazine, found_fluorobenzyl, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product = Chem.MolFromSmiles(product_smiles)

            if product:
                # Check for piperazine scaffold
                piperazine_pattern = Chem.MolFromSmarts("N1CCN(CC1)")
                if product.HasSubstructMatch(piperazine_pattern):
                    found_piperazine = True

                # Check for fluorobenzyl group
                fluorobenzyl_pattern = Chem.MolFromSmarts("c1cc(F)ccc1C")
                if product.HasSubstructMatch(fluorobenzyl_pattern):
                    found_fluorobenzyl = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found a linear synthesis (3+ reactions) with piperazine and fluorobenzyl
    return reaction_count >= 3 and found_piperazine and found_fluorobenzyl
