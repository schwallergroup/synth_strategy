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
    This function detects a strategy involving late-stage functional group modifications,
    particularly focusing on ester hydrolysis in the final steps.
    """
    late_stage_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_hydrolysis

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider reactions at depth 0 or 1 (late stage)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles if r
                    ]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if all(reactant_mols) and product_mol:
                        # Check for ester pattern in reactants
                        ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#6]")
                        has_ester = any(
                            mol.HasSubstructMatch(ester_pattern)
                            for mol in reactant_mols
                        )

                        # Check for carboxylic acid pattern in product
                        acid_pattern = Chem.MolFromSmarts("[#8]-[#6](=[#8])-[#6]")
                        has_acid = product_mol.HasSubstructMatch(acid_pattern)

                        if has_ester and has_acid:
                            print(
                                f"Detected late-stage ester hydrolysis at depth {depth}"
                            )
                            late_stage_hydrolysis = True
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_hydrolysis
