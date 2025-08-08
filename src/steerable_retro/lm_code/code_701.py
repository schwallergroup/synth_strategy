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
    This function detects a functional group interconversion sequence
    from ester to acid to amide throughout the synthesis.
    """
    # Track the sequence of functional groups
    fg_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Define patterns
                ester_pattern = Chem.MolFromSmarts("[C:1](=[O:2])[O][C]")
                acid_pattern = Chem.MolFromSmarts("[C:1](=[O:2])[OH]")
                amide_pattern = Chem.MolFromSmarts("[C:1](=[O:2])[NH]")

                # Check product for functional groups
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol:
                    if product_mol.HasSubstructMatch(ester_pattern):
                        fg_sequence.append(("ester", depth))
                    elif product_mol.HasSubstructMatch(acid_pattern):
                        fg_sequence.append(("acid", depth))
                    elif product_mol.HasSubstructMatch(amide_pattern):
                        fg_sequence.append(("amide", depth))

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort by depth in descending order (from early to late stages)
    fg_sequence.sort(key=lambda x: x[1], reverse=True)

    # Extract just the functional group names
    fg_names = [fg for fg, _ in fg_sequence]

    # Check if the sequence contains ester -> acid -> amide
    has_sequence = False
    if len(fg_names) >= 3:
        for i in range(len(fg_names) - 2):
            if fg_names[i] == "ester" and fg_names[i + 1] == "acid" and fg_names[i + 2] == "amide":
                has_sequence = True
                print(f"Found ester -> acid -> amide sequence: {fg_sequence}")
                break

    return has_sequence
