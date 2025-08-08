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
    This function detects a strategy where a nitrile group is transformed
    into an amide during the synthesis.
    """
    nitrile_to_amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_to_amide_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant has nitrile and product has amide
                reactant_has_nitrile = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
                        if mol.HasSubstructMatch(nitrile_pattern):
                            reactant_has_nitrile = True
                            break

                product_has_amide = False
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
                    if product_mol.HasSubstructMatch(amide_pattern):
                        product_has_amide = True

                if reactant_has_nitrile and product_has_amide:
                    print(f"Nitrile to amide transformation detected at depth {depth}")
                    nitrile_to_amide_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return nitrile_to_amide_found
