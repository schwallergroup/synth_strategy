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
    This function detects a synthetic strategy where a pyrazolone ring is formed
    in a late stage of the synthesis through reaction with hydrazine.
    """
    # Track if we found a pyrazolone formation
    found_pyrazolone_formation = False
    # Track the depth at which it occurs
    formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_pyrazolone_formation, formation_depth

        if node["type"] == "reaction":
            # Check if this is a pyrazolone formation reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if hydrazine is one of the reactants
                hydrazine_pattern = re.compile(r"\[NH2\]\[NH2\]|\[NH2:.*\]\[NH2:.*\]")
                has_hydrazine = any(hydrazine_pattern.search(r) for r in reactants)

                # Check if product contains a pyrazolone structure
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    pyrazolone_pattern = Chem.MolFromSmarts("[#6]1=[#7][#7H][#6](=[#8])[#6]1")
                    if product_mol.HasSubstructMatch(pyrazolone_pattern) and has_hydrazine:
                        found_pyrazolone_formation = True
                        formation_depth = depth
                        print(f"Found pyrazolone formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Consider it late-stage if it occurs at depth 0 or 1
    is_late_stage = formation_depth is not None and formation_depth <= 1

    if found_pyrazolone_formation and is_late_stage:
        print("Detected late-stage pyrazolone formation strategy")
        return True
    return False
