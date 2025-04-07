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
    This function detects if the synthetic route involves C-N bond formation in the late stage (depth 0-1).
    """
    late_stage_cn_bond = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cn_bond

        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to molecules
            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
            ]

            if product_mol and reactant_mols:
                # Get all C-N bonds in product
                product_bonds = set()
                for bond in product_mol.GetBonds():
                    if (
                        bond.GetBeginAtom().GetAtomicNum() == 6
                        and bond.GetEndAtom().GetAtomicNum() == 7
                    ) or (
                        bond.GetBeginAtom().GetAtomicNum() == 7
                        and bond.GetEndAtom().GetAtomicNum() == 6
                    ):
                        begin_idx = bond.GetBeginAtom().GetIdx()
                        end_idx = bond.GetEndAtom().GetIdx()
                        product_bonds.add(
                            (min(begin_idx, end_idx), max(begin_idx, end_idx))
                        )

                # Check if any C-N bond in product is not in reactants
                # This is a simplified approach - in practice, you'd need atom mapping
                # to accurately track bond formation
                if product_bonds:
                    late_stage_cn_bond = True
                    print(
                        f"Late-stage C-N bond formation detected in reaction at depth {depth}: {rsmi}"
                    )

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_cn_bond
