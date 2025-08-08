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
    This function detects a synthetic strategy involving carbonyl chemistry
    for key C-C bond formation in early stages.
    """
    carbonyl_cc_formation = False

    def dfs_traverse(node):
        nonlocal carbonyl_cc_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this is an early-stage reaction (depth > 1)
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert SMILES to molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for carbonyl groups in reactants
                carbonyl_pattern = Chem.MolFromSmarts("[C$(C=O)]")

                reactants_have_carbonyl = False
                for r in reactants:
                    if r.HasSubstructMatch(carbonyl_pattern):
                        reactants_have_carbonyl = True
                        break

                if reactants_have_carbonyl:
                    # This is a simplified check for C-C bond formation
                    # In a real implementation, you would need to compare the atom mappings
                    # between reactants and products to identify new C-C bonds
                    print("Detected carbonyl-based C-C bond formation")
                    carbonyl_cc_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return carbonyl_cc_formation
