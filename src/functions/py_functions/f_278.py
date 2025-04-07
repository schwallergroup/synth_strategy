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
    This function detects if the synthetic route involves a late-stage amide formation.
    """
    late_stage_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check only late-stage reactions (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acyl chloride in reactants
                acyl_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")
                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("C(=O)N")

                acyl_chloride_found = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(acyl_chloride_pattern):
                            acyl_chloride_found = True
                            break
                    except:
                        continue

                if acyl_chloride_found:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(amide_pattern):
                            print("Late-stage amide formation detected at depth", depth)
                            late_stage_amide_formation = True
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_amide_formation
