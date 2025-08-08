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
    This function detects a strategy where the final step (or one of the last steps)
    is an acylation reaction, particularly on a nitrogen atom.
    """
    acylation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal acylation_depths

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acyl chloride pattern in reactants
                acyl_chloride_pattern = Chem.MolFromSmarts("[C](=O)Cl")

                # Check for amide formation in product
                amide_pattern = Chem.MolFromSmarts("[#7][C](=O)")

                acyl_chloride_present = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(acyl_chloride_pattern):
                        acyl_chloride_present = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    acyl_chloride_present
                    and product_mol
                    and product_mol.HasSubstructMatch(amide_pattern)
                ):
                    print(f"Acylation reaction detected at depth {depth}")
                    acylation_depths.append(depth)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if acylation occurs at depth 0 or 1 (late stage)
    return any(depth <= 1 for depth in acylation_depths)
