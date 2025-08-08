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
    Detects incorporation of halogenated aryl groups (like 2,4-dichlorobenzyl)
    via alkylation reactions.
    """
    found_halogenated_aryl_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_halogenated_aryl_incorporation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            if product is None or any(r is None for r in reactants):
                print(f"Warning: Could not parse SMILES at depth {depth}")
                return

            # Check for halogenated aryl group in reactants
            halogenated_aryl_pattern = Chem.MolFromSmarts("[c]([Cl,Br,F,I])")
            alkyl_halide_pattern = Chem.MolFromSmarts("[C]-[Br,Cl,I,F]")

            # Look for dichlorobenzyl specifically
            dichlorobenzyl_pattern = Chem.MolFromSmarts("[#6]-[c]1[c]([Cl])[c][c][c][c]1[Cl]")

            for r in reactants:
                if r.HasSubstructMatch(halogenated_aryl_pattern) and r.HasSubstructMatch(
                    alkyl_halide_pattern
                ):
                    # Check if product now has the halogenated aryl group
                    if product.HasSubstructMatch(halogenated_aryl_pattern):
                        found_halogenated_aryl_incorporation = True
                        print(f"Found halogenated aryl incorporation at depth {depth}")

                        # Check specifically for 2,4-dichlorobenzyl
                        if r.HasSubstructMatch(dichlorobenzyl_pattern):
                            print(f"Specifically found 2,4-dichlorobenzyl incorporation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if found_halogenated_aryl_incorporation:
        print("Found halogenated aryl incorporation strategy")
    else:
        print("Did not find halogenated aryl incorporation")

    return found_halogenated_aryl_incorporation
