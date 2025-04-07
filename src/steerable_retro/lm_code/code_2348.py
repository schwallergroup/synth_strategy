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
    Detects if the synthesis includes an aryl C-N bond formation step.
    """
    aryl_cn_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal aryl_cn_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                reaction_smiles = node["metadata"]["rsmi"]
                reactants_part = reaction_smiles.split(">")[0]
                products_part = reaction_smiles.split(">")[-1]

                # Check for aryl bromide in reactants
                aryl_bromide_pattern = Chem.MolFromSmarts("c[Br]")
                # Check for aniline in reactants
                aniline_pattern = Chem.MolFromSmarts("c[NH2]")
                # Check for aryl amine in product
                aryl_amine_pattern = Chem.MolFromSmarts("c[NH]c")

                reactants_mol = Chem.MolFromSmiles(reactants_part)
                product_mol = Chem.MolFromSmiles(products_part)

                if (
                    reactants_mol
                    and product_mol
                    and reactants_mol.HasSubstructMatch(aryl_bromide_pattern)
                    and reactants_mol.HasSubstructMatch(aniline_pattern)
                    and product_mol.HasSubstructMatch(aryl_amine_pattern)
                ):
                    aryl_cn_formation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Aryl C-N bond formation: {aryl_cn_formation}")

    return aryl_cn_formation
