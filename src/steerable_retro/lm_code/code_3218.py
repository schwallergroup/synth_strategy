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
    Detects if the synthesis route involves a late-stage vinyl coupling reaction
    (introduction of a C=C bond in the final or penultimate step).
    """
    vinyl_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal vinyl_coupling_detected

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if vinyl group is being introduced
                vinyl_pattern = Chem.MolFromSmarts("[C]=[C]")

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(vinyl_pattern):
                    # Check if vinyl group was not present in all reactants
                    vinyl_in_reactants = all(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(vinyl_pattern)
                        for r in reactants
                        if Chem.MolFromSmiles(r)
                    )

                    if not vinyl_in_reactants:
                        print(f"Late-stage vinyl coupling detected at depth {depth}")
                        vinyl_coupling_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return vinyl_coupling_detected
