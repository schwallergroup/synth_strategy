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
    This function detects if the synthesis route employs an ester → acid → amide
    functional group transformation sequence.
    """
    # Track the sequence of functional group transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for functional groups
            ester_pattern = Chem.MolFromSmarts("[C:1](=[O:2])[O:3][C:4]")
            acid_pattern = Chem.MolFromSmarts("[C:1](=[O:2])[OH:3]")
            amide_pattern = Chem.MolFromSmarts("[C:1](=[O:2])[N:3]")

            # Check reactants and products
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol:
                # Ester to acid transformation
                if any(
                    mol.HasSubstructMatch(ester_pattern) for mol in reactant_mols if mol
                ) and product_mol.HasSubstructMatch(acid_pattern):
                    transformations.append(("ester_to_acid", depth))
                    print(f"Found ester to acid transformation at depth {depth}")

                # Acid to amide transformation
                if any(
                    mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols if mol
                ) and product_mol.HasSubstructMatch(amide_pattern):
                    transformations.append(("acid_to_amide", depth))
                    print(f"Found acid to amide transformation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if the sequence ester → acid → amide exists in the correct order
    has_ester_to_acid = any(t[0] == "ester_to_acid" for t in transformations)
    has_acid_to_amide = any(t[0] == "acid_to_amide" for t in transformations)

    # Check if they occur in the correct sequence (ester_to_acid should have higher depth than acid_to_amide)
    correct_sequence = False
    if has_ester_to_acid and has_acid_to_amide:
        ester_to_acid_depth = next(t[1] for t in transformations if t[0] == "ester_to_acid")
        acid_to_amide_depth = next(t[1] for t in transformations if t[0] == "acid_to_amide")
        correct_sequence = ester_to_acid_depth > acid_to_amide_depth

    print(f"Transformations: {transformations}, Has correct sequence: {correct_sequence}")
    return correct_sequence
