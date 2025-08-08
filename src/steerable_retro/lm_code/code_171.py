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
    Detects if the synthetic route involves incorporation of a morpholine group
    in a late-stage reaction.
    """
    morpholine_late_stage = False

    def dfs_traverse(node):
        nonlocal morpholine_late_stage

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            depth = node.get("metadata", {}).get("depth", -1)

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")

                # Check for morpholine pattern in reactants
                morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(morpholine_pattern):
                        # Check if product also has morpholine
                        product_mol = Chem.MolFromSmiles(product_part)
                        if product_mol and product_mol.HasSubstructMatch(morpholine_pattern):
                            print("Found late-stage morpholine incorporation")
                            morpholine_late_stage = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return morpholine_late_stage
