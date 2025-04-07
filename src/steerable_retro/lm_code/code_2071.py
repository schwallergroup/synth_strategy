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
    This function detects if the synthesis involves late-stage amide formation.
    """
    late_amide_formation = False
    low_depth_threshold = 1  # Depth <= 1 is considered late stage

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation

        if node["type"] == "reaction" and depth <= low_depth_threshold:
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check for carboxylic acid in reactants
                acid_pattern = Chem.MolFromSmarts("[#6](=[O])-[O;H1]")

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[#6](=[O])-[#7]")

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    acid_found = False
                    amine_found = False

                    for r_mol in reactant_mols:
                        if r_mol:
                            if r_mol.HasSubstructMatch(acid_pattern):
                                acid_found = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                amine_found = True

                    if acid_found and amine_found:
                        print(f"Found late-stage amide formation at depth {depth}")
                        late_amide_formation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_amide_formation
