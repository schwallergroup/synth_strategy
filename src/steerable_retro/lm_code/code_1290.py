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
    This function detects if the synthesis involves a late-stage amide formation
    (in the last 40% of steps)
    """
    # Track amide formation reactions and their depths
    amide_formation_depths = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an amide formation reaction
            product_mol = Chem.MolFromSmiles(product)
            has_amide_product = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[#6]C(=O)[N][#6]")
            )

            # Check if reactants have carboxylic acid and amine
            has_acid = False
            has_amine = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]C(=O)[OH]")):
                        has_acid = True
                    if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2][#6]")):
                        has_amine = True

            # If this looks like amide formation, record the depth
            if has_amide_product and (has_acid or has_amine):
                amide_formation_depths.append(depth)
                print(f"Found amide formation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if amide formation occurs in the last 40% of steps
    if amide_formation_depths and max_depth > 0:
        # Lower depths are later in the synthesis (closer to final product)
        late_stage_threshold = 0.4 * max_depth
        is_late_stage = min(amide_formation_depths) <= late_stage_threshold
        print(f"Late stage amide formation detected: {is_late_stage}")
        return is_late_stage

    print("No amide formation detected")
    return False
