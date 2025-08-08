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
    Detects a sequence of functional group interconversions:
    ester → alcohol → halide for activation.
    """
    # Initialize tracking variables
    ester_to_alcohol = False
    alcohol_to_halide = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_to_alcohol, alcohol_to_halide

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                # Check for ester to alcohol reduction
                if product:
                    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
                    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")

                    has_ester_reactant = any(
                        r and ester_pattern and r.HasSubstructMatch(ester_pattern)
                        for r in reactants
                    )
                    has_alcohol_product = (
                        product and alcohol_pattern and product.HasSubstructMatch(alcohol_pattern)
                    )

                    if has_ester_reactant and has_alcohol_product:
                        ester_to_alcohol = True
                        print("Detected ester to alcohol reduction")

                # Check for alcohol to halide conversion
                if product:
                    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
                    halide_pattern = Chem.MolFromSmarts("[CX4][Br,Cl,I,F]")

                    has_alcohol_reactant = any(
                        r and alcohol_pattern and r.HasSubstructMatch(alcohol_pattern)
                        for r in reactants
                    )
                    has_halide_product = (
                        product and halide_pattern and product.HasSubstructMatch(halide_pattern)
                    )

                    if has_alcohol_reactant and has_halide_product:
                        alcohol_to_halide = True
                        print("Detected alcohol to halide conversion")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the complete sequence is present
    sequence_present = ester_to_alcohol and alcohol_to_halide

    if sequence_present:
        print(
            "Detected complete functional group interconversion sequence: ester → alcohol → halide"
        )

    return sequence_present
