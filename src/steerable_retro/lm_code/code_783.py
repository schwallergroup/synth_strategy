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
    Detects nucleophilic substitution where morpholine displaces a leaving group.
    """
    morpholine_substitution_detected = False

    def dfs_traverse(node):
        nonlocal morpholine_substitution_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for morpholine pattern in reactants
            morpholine_pattern = Chem.MolFromSmarts("[NH]1CCOCC1")
            has_morpholine = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(morpholine_pattern):
                    has_morpholine = True
                    break

            # Check for leaving group pattern in reactants (mesylate, tosylate, halide)
            leaving_group_pattern = Chem.MolFromSmarts("[C]-[O,S,F,Cl,Br,I]")
            has_leaving_group = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(leaving_group_pattern):
                    has_leaving_group = True
                    break

            # Check if product has C-N bond to morpholine
            morpholine_product_pattern = Chem.MolFromSmarts("[C]-[N]1CCOCC1")
            product_mol = Chem.MolFromSmiles(product)

            if (
                has_morpholine
                and has_leaving_group
                and product_mol
                and product_mol.HasSubstructMatch(morpholine_product_pattern)
            ):
                print("Nucleophilic substitution with morpholine detected")
                morpholine_substitution_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return morpholine_substitution_detected
