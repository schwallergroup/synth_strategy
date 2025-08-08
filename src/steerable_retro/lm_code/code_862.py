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
    This function detects C-N bond formations occurring in the second half of the synthesis.
    """
    late_stage_cn_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cn_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # In retrosynthetic direction, low depth means late stage
            if depth <= 3:  # Second half of synthesis (assuming ~7 steps total)
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C-N bond formation
                for reactant in reactants:
                    if not reactant:
                        continue

                    r_mol = Chem.MolFromSmiles(reactant)
                    p_mol = Chem.MolFromSmiles(product)

                    if r_mol is None or p_mol is None:
                        continue

                    # Count C-N bonds in reactant and product
                    cn_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                    r_cn_count = len(r_mol.GetSubstructMatches(cn_pattern))
                    p_cn_count = len(p_mol.GetSubstructMatches(cn_pattern))

                    if p_cn_count > r_cn_count:
                        late_stage_cn_formation = True
                        print(
                            f"Found late-stage C-N bond formation at depth {depth}, reaction: {rsmi}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_stage_cn_formation
