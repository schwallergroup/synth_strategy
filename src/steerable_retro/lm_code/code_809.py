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
    This function detects a late-stage functionalization strategy, particularly
    with trifluoromethyl-containing groups.
    """
    late_stage_functionalization = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_functionalization

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for trifluoromethyl group in product
            trifluoromethyl_pattern = Chem.MolFromSmarts("[#6]([#9])([#9])[#9]")

            # Check if this is a late-stage reaction (depth <= 2)
            if depth <= 2:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(trifluoromethyl_pattern):
                    # Check if any reactant doesn't have trifluoromethyl
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and not mol.HasSubstructMatch(trifluoromethyl_pattern):
                            late_stage_functionalization = True
                            print(
                                f"Detected late-stage functionalization with trifluoromethyl group: {rsmi}"
                            )
                            break

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage functionalization strategy detected: {late_stage_functionalization}")
    return late_stage_functionalization
