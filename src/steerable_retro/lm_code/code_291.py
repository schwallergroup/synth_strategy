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
    This function detects late-stage sulfonamide formation from an amine,
    typically in the final steps of the synthesis.
    """
    found_late_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_sulfonamide

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonamide formation
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol:
                    # Check if product contains sulfonamide
                    sulfonamide_pattern = Chem.MolFromSmarts("[#7]S(=O)(=O)[#6]")

                    # Check if any reactant contains primary amine
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    amine_pattern = Chem.MolFromSmarts("[#7;H2]")

                    reactants_have_amine = False
                    for mol in reactant_mols:
                        if mol and mol.HasSubstructMatch(amine_pattern):
                            reactants_have_amine = True
                            break

                    if product_mol.HasSubstructMatch(sulfonamide_pattern) and reactants_have_amine:
                        found_late_sulfonamide = True
                        print(f"Found late-stage sulfonamide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_late_sulfonamide
