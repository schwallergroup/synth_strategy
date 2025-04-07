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
    Detects an ester hydrolysis in the final step of the synthesis.
    """
    final_step_is_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_hydrolysis

        if node["type"] == "reaction" and depth == 0:  # Final step (depth 0)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if (
                        all(mol is not None for mol in reactant_mols)
                        and product_mol is not None
                    ):
                        # Check for ester hydrolysis
                        ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")
                        acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")

                        has_ester = any(
                            mol.HasSubstructMatch(ester_pattern)
                            for mol in reactant_mols
                        )
                        has_acid = product_mol.HasSubstructMatch(acid_pattern)

                        if has_ester and has_acid:
                            final_step_is_hydrolysis = True
                            print("Ester hydrolysis detected in final step")
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return final_step_is_hydrolysis
