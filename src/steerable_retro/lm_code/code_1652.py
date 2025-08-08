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
    Detects if the synthesis follows a linear strategy (each reaction has only one product).

    In a linear synthesis:
    1. Each reaction produces exactly one product molecule
    2. Each intermediate molecule is used in exactly one reaction
    3. The final target molecule is not used as a reactant
    """
    # Track how many times each molecule is used as a reactant
    molecule_usage = {}
    reaction_count = 0

    # First pass: count molecule usage
    def count_usage(node, parent_type=None):
        nonlocal reaction_count

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Initialize if not seen before
            if mol_smiles not in molecule_usage:
                molecule_usage[mol_smiles] = {"as_reactant": 0, "as_product": 0}

            # If parent was a reaction, this molecule is a reactant in retrosynthesis
            if parent_type == "reaction":
                molecule_usage[mol_smiles]["as_reactant"] += 1

        elif node["type"] == "reaction":
            reaction_count += 1
            # Find the parent molecule (product in retrosynthesis)
            for child in node.get("children", []):
                if child["type"] == "mol":
                    mol_smiles = child["smiles"]
                    if mol_smiles not in molecule_usage:
                        molecule_usage[mol_smiles] = {"as_reactant": 0, "as_product": 0}
                    molecule_usage[mol_smiles]["as_product"] += 1

        # Traverse children
        for child in node.get("children", []):
            count_usage(child, node["type"])

    # Second pass: check if synthesis is linear
    def check_linearity():
        # No reactions means not a synthesis
        if reaction_count == 0:
            print("No reactions found in the route")
            return False

        # Check if each molecule follows linear pattern
        for smiles, usage in molecule_usage.items():
            # Starting materials should not be products
            if usage["as_reactant"] == 0 and usage["as_product"] == 0:
                continue  # Skip in-stock molecules with no reactions

            # Target molecule: used as product once, never as reactant
            if usage["as_product"] == 0 and usage["as_reactant"] == 1:
                continue

            # Intermediate molecules: used once as product, once as reactant
            if usage["as_product"] == 1 and usage["as_reactant"] == 1:
                continue

            # Any other pattern means non-linear synthesis
            print(f"Non-linear pattern found for molecule: {smiles}, usage: {usage}")
            return False

        return True

    # Execute the analysis
    count_usage(route)
    result = check_linearity()

    print(f"Linear synthesis strategy: {result}")
    return result
