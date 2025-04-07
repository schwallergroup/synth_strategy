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
    This function detects a synthetic strategy involving amide formation
    using an acyl chloride intermediate.
    """
    has_acyl_chloride_amidation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_acyl_chloride_amidation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acyl chloride pattern in reactants
                acyl_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")

                # Check for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                acyl_chloride_present = False
                for r in reactants:
                    try:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(acyl_chloride_pattern):
                            acyl_chloride_present = True
                            break
                    except:
                        continue

                if acyl_chloride_present:
                    try:
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol and p_mol.HasSubstructMatch(amide_pattern):
                            print(f"Detected amide formation via acyl chloride at depth {depth}")
                            has_acyl_chloride_amidation = True
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Strategy amide_formation_via_acyl_chloride: {has_acyl_chloride_amidation}")
    return has_acyl_chloride_amidation
