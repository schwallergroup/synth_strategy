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
    Detects if the synthesis involves sulfonamide formation in the late stage.
    Late stage is defined as low depth in the retrosynthetic tree.
    """
    sulfonamide_formed_late = False
    late_stage_threshold = 1  # Consider depth <= 1 as late stage

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formed_late

        if node["type"] == "reaction" and depth <= late_stage_threshold:
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this reaction forms a sulfonamide
                product_mol = Chem.MolFromSmiles(product)
                sulfonamide_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[N]")

                if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                    # Check if reactants include sulfonyl chloride
                    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[Cl]")
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            sulfonyl_chloride_pattern
                        ):
                            print(f"Late-stage sulfonamide formation detected at depth {depth}")
                            sulfonamide_formed_late = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return sulfonamide_formed_late
