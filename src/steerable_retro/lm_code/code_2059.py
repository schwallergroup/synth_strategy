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
    Detects sulfonamide formation in the synthesis route.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonyl chloride in reactants
                sulfonyl_chloride_present = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[S](=O)(=O)[Cl]")

                    if (
                        reactant_mol
                        and sulfonyl_chloride_pattern
                        and reactant_mol.HasSubstructMatch(sulfonyl_chloride_pattern)
                    ):
                        sulfonyl_chloride_present = True
                        break

                # Check for sulfonamide in product
                if sulfonyl_chloride_present:
                    product_mol = Chem.MolFromSmiles(product)
                    sulfonamide_pattern = Chem.MolFromSmarts("[S](=O)(=O)[N]")

                    if (
                        product_mol
                        and sulfonamide_pattern
                        and product_mol.HasSubstructMatch(sulfonamide_pattern)
                    ):
                        has_sulfonamide_formation = True
                        print(
                            f"Detected sulfonamide formation at depth {node.get('depth', 'unknown')}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_sulfonamide_formation
