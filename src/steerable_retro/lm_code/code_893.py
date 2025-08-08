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
    Detects a synthetic strategy involving late-stage halogenation (iodination, bromination).
    Late stage means in the first third of the synthesis (lower depth values).
    """
    late_halogenation_found = False
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal late_halogenation_found, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(r for r in reactants):
                # Check for halogenation (C-H to C-X transformation)
                halogen_pattern = Chem.MolFromSmarts("[c]-[I,Br]")

                product_has_halogen = product.HasSubstructMatch(halogen_pattern)
                reactants_have_halogen = any(
                    r.HasSubstructMatch(halogen_pattern) for r in reactants
                )

                if (
                    product_has_halogen and not reactants_have_halogen and depth <= 2
                ):  # Late stage (depth ≤ 2)
                    late_halogenation_found = True
                    print(f"Late-stage halogenation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage halogenation strategy: {late_halogenation_found}")
    return late_halogenation_found
