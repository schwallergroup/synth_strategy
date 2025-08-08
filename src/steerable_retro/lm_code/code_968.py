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
    This function detects if the synthetic route includes a late-stage introduction
    of morpholine via nucleophilic substitution.
    """
    morpholine_introduced = False

    def dfs_traverse(node):
        nonlocal morpholine_introduced

        # Check only reactions at depth 0 or 1 (late-stage)
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for morpholine in reactants
            morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#6]([#6][#6]1)[#8]")
            # Check for chloromethyl in reactants
            chloromethyl_pattern = Chem.MolFromSmarts("[c][CH2][Cl]")
            # Check for morpholine attachment in product
            morpholine_attached_pattern = Chem.MolFromSmarts(
                "[c][CH2][#7]1[#6][#6][#6]([#6][#6]1)[#8]"
            )

            morpholine_present = False
            chloromethyl_present = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(morpholine_pattern):
                        morpholine_present = True
                    if mol and mol.HasSubstructMatch(chloromethyl_pattern):
                        chloromethyl_present = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                morpholine_attached = product_mol and product_mol.HasSubstructMatch(
                    morpholine_attached_pattern
                )
            except:
                morpholine_attached = False

            if morpholine_present and chloromethyl_present and morpholine_attached:
                print("Late-stage morpholine introduction detected")
                morpholine_introduced = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return morpholine_introduced
