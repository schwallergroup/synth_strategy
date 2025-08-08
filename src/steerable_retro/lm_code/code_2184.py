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
    This function detects if N-benzylation of indole occurs in the synthesis.
    """
    n_benzylation_detected = False

    def dfs_traverse(node):
        nonlocal n_benzylation_detected
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for indole N-H in reactants
                indole_nh_pattern = Chem.MolFromSmarts("[#6]1[#6][#7H][#6]2[#6][#6][#6][#6][#6]12")
                # Check for N-benzylated indole in product
                n_benzyl_indole_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#7]([#6][#6]2[#6][#6][#6][#6][#6]2)[#6]3[#6][#6][#6][#6][#6]13"
                )

                reactant_has_indole_nh = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(indole_nh_pattern):
                        reactant_has_indole_nh = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                product_has_n_benzyl = product_mol and product_mol.HasSubstructMatch(
                    n_benzyl_indole_pattern
                )

                if reactant_has_indole_nh and product_has_n_benzyl:
                    n_benzylation_detected = True
                    print(f"N-benzylation detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"N-benzylation of indole detected: {n_benzylation_detected}")
    return n_benzylation_detected
