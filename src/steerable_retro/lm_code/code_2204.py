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
    This function detects if the synthesis involves a sequence of ester → acid → amide conversions.
    """
    # Track the sequence of functional groups in reverse synthetic direction
    functional_group_sequence = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol is not None:
                    # Check for functional groups
                    ester_pattern = Chem.MolFromSmarts("[#6][#8][#6](=[#8])")
                    acid_pattern = Chem.MolFromSmarts("[#8H][#6](=[#8])")
                    amide_pattern = Chem.MolFromSmarts("[#6][#6](=[#8])[#7]")

                    if product_mol.HasSubstructMatch(ester_pattern):
                        functional_group_sequence.append("ester")
                    elif product_mol.HasSubstructMatch(acid_pattern):
                        functional_group_sequence.append("acid")
                    elif product_mol.HasSubstructMatch(amide_pattern):
                        functional_group_sequence.append("amide")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if the sequence contains ester → acid → amide (in reverse order since we're traversing retrosynthetically)
    has_sequence = False
    for i in range(len(functional_group_sequence) - 2):
        if (
            functional_group_sequence[i] == "amide"
            and functional_group_sequence[i + 1] == "acid"
            and functional_group_sequence[i + 2] == "ester"
        ):
            has_sequence = True
            print("Detected ester → acid → amide conversion sequence")
            break

    return has_sequence
