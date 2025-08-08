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
    Detects Suzuki coupling (aryl halide + boronate) in the synthesis route.
    """
    found_suzuki = False

    def dfs_traverse(node):
        nonlocal found_suzuki

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                # Check for boronate pattern
                boronate_pattern = Chem.MolFromSmarts("[#6]B([O])[O]")
                has_boronate = any(
                    mol.GetSubstructMatches(boronate_pattern) for mol in reactants if mol
                )

                # Check for aryl halide pattern
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I]")
                has_aryl_halide = any(
                    mol.GetSubstructMatches(aryl_halide_pattern) for mol in reactants if mol
                )

                if has_boronate and has_aryl_halide:
                    print("Found Suzuki coupling pattern")
                    found_suzuki = True
            except:
                print("Error processing reaction SMILES for Suzuki detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_suzuki
