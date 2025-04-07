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
    This function detects if the synthetic route involves late-stage SNAr to form diaryl ether.
    """
    diaryl_ether_formed = False
    formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal diaryl_ether_formed, formation_depth

        if node["type"] == "reaction":
            # Extract product and reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for diaryl ether formation
            diaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                # Check if reactants have phenol and aryl halide
                phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I,F]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                has_phenol = any(
                    mol and mol.HasSubstructMatch(phenol_pattern)
                    for mol in reactant_mols
                )
                has_aryl_halide = any(
                    mol and mol.HasSubstructMatch(aryl_halide_pattern)
                    for mol in reactant_mols
                )

                if has_phenol and has_aryl_halide:
                    diaryl_ether_formed = True
                    formation_depth = depth
                    print(f"Diaryl ether formation via SNAr detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Late stage is defined as depth <= 1 (closer to final product)
    return diaryl_ether_formed and formation_depth <= 1
