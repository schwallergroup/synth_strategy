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
    This function detects the use of trifluoroacetyl as a protecting group
    that is maintained throughout multiple steps of the synthesis.
    """
    # Track reactions containing trifluoroacetyl group
    reactions_with_tfa = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product_smiles = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol is not None:
                    # Check for trifluoroacetyl group
                    tfa_pattern = Chem.MolFromSmarts("[#6](=[#8])[#6]([#9])([#9])[#9]")
                    if product_mol.HasSubstructMatch(tfa_pattern):
                        reactions_with_tfa.append(depth)
                        print(f"Trifluoroacetyl group detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if trifluoroacetyl group is present in multiple consecutive reactions
    if len(reactions_with_tfa) >= 3:
        print(
            f"Trifluoroacetyl protection strategy detected across {len(reactions_with_tfa)} reactions"
        )
        return True
    return False
