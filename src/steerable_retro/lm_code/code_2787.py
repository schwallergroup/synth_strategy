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
    This function detects a convergent synthesis with late-stage amide coupling.
    It looks for a reaction at low depth (late stage) that forms an amide bond
    and combines multiple fragments.
    """
    late_stage_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage reaction (depth 0 or 1)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if we have multiple reactants (convergent)
                if len(reactants_smiles) >= 2:
                    # Check for amide formation
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Count amide bonds in reactants and product
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                    reactant_amide_count = sum(
                        [
                            len(mol.GetSubstructMatches(amide_pattern)) if mol else 0
                            for mol in reactant_mols
                        ]
                    )
                    product_amide_count = (
                        len(product_mol.GetSubstructMatches(amide_pattern)) if product_mol else 0
                    )

                    # If product has more amide bonds than reactants combined, amide formation occurred
                    if product_amide_count > reactant_amide_count:
                        print(f"Found late-stage amide coupling at depth {depth}")
                        late_stage_amide_coupling = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_amide_coupling
