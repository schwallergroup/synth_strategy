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
    This function detects if the synthesis involves early heterocycle formation
    followed by late-stage functionalization.
    """
    heterocycle_depth = None
    functionalization_depths = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_depth, functionalization_depths, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                # Check for heterocycle formation (pyrazole)
                pyrazole_patt = Chem.MolFromSmarts("c1nn[c]c1")
                if product_mol.HasSubstructMatch(pyrazole_patt):
                    # Check if reactants don't have pyrazole
                    pyrazole_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            pyrazole_patt
                        ):
                            pyrazole_in_reactants = True

                    if not pyrazole_in_reactants:
                        heterocycle_depth = depth
                        print(f"Detected heterocycle formation at depth {depth}")

                # Check for late functionalizations
                # 1. Ether formation
                ether_patt = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
                # 2. Cyanation
                nitrile_patt = Chem.MolFromSmarts("C#N")
                # 3. Formylation
                aldehyde_patt = Chem.MolFromSmarts("[CX3H1](=O)")

                for patt, name in [
                    (ether_patt, "ether"),
                    (nitrile_patt, "nitrile"),
                    (aldehyde_patt, "aldehyde"),
                ]:
                    if product_mol.HasSubstructMatch(patt):
                        # Check if reactants don't have the pattern
                        pattern_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(patt):
                                pattern_in_reactants = True

                        if not pattern_in_reactants:
                            functionalization_depths.append(depth)
                            print(f"Detected {name} functionalization at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if heterocycle formation is early (higher depth) and functionalizations are late (lower depth)
    if heterocycle_depth is not None and functionalization_depths:
        early_heterocycle = heterocycle_depth > max_depth / 2
        late_functionalizations = any(
            depth < max_depth / 2 for depth in functionalization_depths
        )

        if early_heterocycle and late_functionalizations:
            print(
                f"Confirmed early heterocycle formation (depth {heterocycle_depth}) with late functionalizations"
            )
            return True

    return False
