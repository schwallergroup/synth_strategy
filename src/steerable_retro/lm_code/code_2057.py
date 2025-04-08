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
    This function detects a strategy involving late-stage introduction of a sulfonyl group.
    """
    late_stage_sulfonyl = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_sulfonyl

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and all(r is not None for r in reactant_mols):
                # Check for sulfonyl group introduction
                sulfonyl_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])")
                if product_mol.HasSubstructMatch(sulfonyl_pattern):
                    # Check if sulfonyl group is coming from one of the reactants
                    sulfonyl_reactant = None
                    for r in reactant_mols:
                        if r.HasSubstructMatch(sulfonyl_pattern):
                            sulfonyl_reactant = r
                            break

                    if sulfonyl_reactant:
                        # Check if the sulfonyl-containing reactant is being incorporated
                        # This is a simplified check - would need more sophisticated matching in practice
                        print(f"Late-stage sulfonyl introduction detected at depth {depth}")
                        late_stage_sulfonyl = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Strategy detected: {late_stage_sulfonyl}")
    return late_stage_sulfonyl
