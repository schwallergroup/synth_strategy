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
    This function detects a synthetic strategy involving a nitrile group,
    particularly looking for late-stage nitrile modifications.
    """
    # Track if we found the pattern
    found_pattern = False
    # Track nitrile-related transformations
    nitrile_present = False
    nitrile_modification = False

    def dfs_traverse(node):
        nonlocal found_pattern, nitrile_present, nitrile_modification

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrile in product
            prod_mol = Chem.MolFromSmiles(product)
            if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]#[#7]")):
                nitrile_present = True
                print(f"Found nitrile in product: {product}")

                # Check for nitrile formation or modification
                if not any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[#6]#[#7]"))
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                ):
                    nitrile_modification = True
                    print(f"Found nitrile formation/modification: {rsmi}")

                    # Check specifically for amidoxime to nitrile conversion
                    if any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(
                            Chem.MolFromSmarts("[#6](=[#7]-[#8])-[#7]")
                        )
                        for r in reactants
                        if Chem.MolFromSmiles(r)
                    ):
                        print(f"Found amidoxime to nitrile conversion: {rsmi}")
                        found_pattern = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If we didn't find the specific amidoxime pattern but found nitrile presence
    if nitrile_present and not found_pattern:
        found_pattern = nitrile_present

    return found_pattern
