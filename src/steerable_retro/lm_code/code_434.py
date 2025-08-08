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
    This function detects if the synthesis involves late-stage nitrile formation
    from an oxime intermediate.
    """
    # Initialize tracking variables
    has_late_stage_nitrile_formation = False

    def dfs_traverse(node):
        nonlocal has_late_stage_nitrile_formation

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 0 or 1)
            depth = node.get("metadata", {}).get("depth", -1)
            if depth <= 1:
                # Check for oxime to nitrile transformation
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    oxime_pattern = Chem.MolFromSmarts("[C]=[N][O]")
                    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")

                    if reactant_mol.HasSubstructMatch(
                        oxime_pattern
                    ) and product_mol.HasSubstructMatch(nitrile_pattern):
                        has_late_stage_nitrile_formation = True
                        print(f"Found late-stage nitrile formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage nitrile formation strategy detected: {has_late_stage_nitrile_formation}")
    return has_late_stage_nitrile_formation
