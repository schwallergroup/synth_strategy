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
    Detects if the synthetic route includes a deprotection as one of the final steps.
    """
    deprotection_depths = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal deprotection_depths, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                # Check for common protecting groups
                boc_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#7]")
                methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6;H3]")

                # Boc deprotection
                if reactant_mol.HasSubstructMatch(
                    boc_pattern
                ) and not product_mol.HasSubstructMatch(boc_pattern):
                    deprotection_depths.append(depth)
                    print(f"Deprotection detected at depth {depth}")

                # O-demethylation
                if reactant_mol.HasSubstructMatch(
                    methoxy_pattern
                ) and not product_mol.HasSubstructMatch(methoxy_pattern):
                    deprotection_depths.append(depth)
                    print(f"Deprotection detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if any deprotection occurs in the last third of the synthesis
    if deprotection_depths and max_depth > 0:
        for depth in deprotection_depths:
            if depth <= max_depth / 3:  # Lower depth values are later in the synthesis
                return True

    return False
