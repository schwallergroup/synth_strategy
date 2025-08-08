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
    This function detects if the synthesis involves hydroxylamine addition to a ketone.
    """
    hydroxylamine_addition_found = False

    def dfs_traverse(node):
        nonlocal hydroxylamine_addition_found

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ketone pattern
            ketone_pattern = Chem.MolFromSmarts("C=O")

            # Check for hydroxylamine pattern
            hydroxylamine_pattern = Chem.MolFromSmarts("[OH][N]")

            # Check for hydroxamic acid pattern
            hydroxamic_acid_pattern = Chem.MolFromSmarts("[OH][N]=C")

            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and product_mol.HasSubstructMatch(hydroxamic_acid_pattern):
                # Check if reactants contain ketone and hydroxylamine
                has_ketone = False
                has_hydroxylamine = False

                for r in reactants:
                    if r:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            if r_mol.HasSubstructMatch(ketone_pattern):
                                has_ketone = True
                            if (
                                r_mol.HasSubstructMatch(hydroxylamine_pattern)
                                or r == "[OH:1][NH2:2]"
                            ):
                                has_hydroxylamine = True

                if has_ketone and has_hydroxylamine:
                    print("Found hydroxylamine addition to ketone")
                    hydroxylamine_addition_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return hydroxylamine_addition_found
