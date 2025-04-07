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
    Detects if the synthetic route involves amide formation in the final step.
    """
    late_amide_formation = False

    def dfs_traverse(node):
        nonlocal late_amide_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this is the final step (depth 0)
            if node.get("depth", 0) == 0:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation
                has_carboxylic_acid = False
                has_amine = False

                for reactant in reactants:
                    if reactant:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]C=O")):
                                has_carboxylic_acid = True
                            if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]")):
                                has_amine = True

                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]C=O")):
                    if has_carboxylic_acid or has_amine:
                        print("Late-stage amide formation detected")
                        late_amide_formation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return late_amide_formation
