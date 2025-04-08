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
    This function detects late-stage amide coupling (amine + carboxylic acid).
    """
    has_late_stage_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage = low depth
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine and carboxylic acid in reactants
                amine_pattern = Chem.MolFromSmarts("c[NH2]")
                carboxylic_pattern = Chem.MolFromSmarts("[#6]C(=O)[OH]")

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")

                has_amine = False
                has_carboxylic = False
                has_amide = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                            if mol.HasSubstructMatch(carboxylic_pattern):
                                has_carboxylic = True
                    except:
                        continue

                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(amide_pattern):
                        has_amide = True
                except:
                    pass

                if has_amine and has_carboxylic and has_amide:
                    print("Detected late-stage amide coupling")
                    has_late_stage_amide_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_amide_coupling
