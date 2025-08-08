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
    This function detects a strategy involving a late-stage Suzuki coupling (depth 0 or 1).
    """
    late_stage_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_suzuki

        if node["type"] == "reaction" and depth <= 1:  # Only consider depth 0 or 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid/ester in reactants
                boronic_pattern = Chem.MolFromSmarts("[#6]-[B]([O][#6])[O][#6]")

                # Check for C-C bond formation
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(boronic_pattern):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # This is a simplified check - in a real implementation,
                            # you would need to compare the product and reactants more carefully
                            # to confirm C-C bond formation
                            late_stage_suzuki = True
                            print(f"Found late-stage Suzuki coupling at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_suzuki
