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
    Detects if the synthetic route employs a transformation of a primary amide
    to a nitrile group.
    """
    found_amide_to_nitrile = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_to_nitrile

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide to nitrile transformation
                primary_amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH2]")
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")

                has_primary_amide = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(primary_amide_pattern):
                            has_primary_amide = True
                            break
                    except:
                        continue

                has_nitrile = False
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(nitrile_pattern):
                        has_nitrile = True
                except:
                    pass

                if has_primary_amide and has_nitrile:
                    found_amide_to_nitrile = True
                    print(f"Found amide to nitrile transformation at depth {depth}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_amide_to_nitrile
