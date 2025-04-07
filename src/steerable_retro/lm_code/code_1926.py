#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects a sequence of ester-acid interconversions in the synthesis route.
    """
    ester_to_acid_count = 0
    acid_to_ester_count = 0

    def dfs_traverse(node):
        nonlocal ester_to_acid_count, acid_to_ester_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester hydrolysis (ester to acid)
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
                acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

                reactant_has_ester = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(ester_pattern):
                        reactant_has_ester = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_ester
                    and product_mol
                    and product_mol.HasSubstructMatch(acid_pattern)
                ):
                    ester_to_acid_count += 1
                    print(f"Detected ester hydrolysis, count: {ester_to_acid_count}")

                # Check for esterification (acid to ester)
                reactant_has_acid = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(acid_pattern):
                        reactant_has_acid = True
                        break

                if (
                    reactant_has_acid
                    and product_mol
                    and product_mol.HasSubstructMatch(ester_pattern)
                ):
                    acid_to_ester_count += 1
                    print(f"Detected esterification, count: {acid_to_ester_count}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if we have at least one of each conversion type
    return ester_to_acid_count >= 1 and acid_to_ester_count >= 1
