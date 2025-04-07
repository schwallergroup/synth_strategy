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
    Detects if the synthesis route involves a late-stage C-N bond formation
    (in the final or penultimate step).
    """
    late_stage_cn_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cn_formation

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Check for C-N bond formation
                    c_n_bonds_product = set()
                    for bond in product.GetBonds():
                        atoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
                        if (
                            bond.GetBeginAtom().GetAtomicNum() == 6
                            and bond.GetEndAtom().GetAtomicNum() == 7
                        ) or (
                            bond.GetBeginAtom().GetAtomicNum() == 7
                            and bond.GetEndAtom().GetAtomicNum() == 6
                        ):
                            idx1, idx2 = sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
                            c_n_bonds_product.add((idx1, idx2))

                    # Check if these C-N bonds exist in reactants
                    c_n_bonds_reactants = set()
                    for r in reactants:
                        for bond in r.GetBonds():
                            atoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
                            if (
                                bond.GetBeginAtom().GetAtomicNum() == 6
                                and bond.GetEndAtom().GetAtomicNum() == 7
                            ) or (
                                bond.GetBeginAtom().GetAtomicNum() == 7
                                and bond.GetEndAtom().GetAtomicNum() == 6
                            ):
                                idx1, idx2 = sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
                                c_n_bonds_reactants.add((idx1, idx2))

                    # If there are more C-N bonds in product than in reactants, a C-N bond was formed
                    if len(c_n_bonds_product) > len(c_n_bonds_reactants):
                        print(f"Late-stage C-N bond formation detected at depth {depth}")
                        late_stage_cn_formation = True
            except:
                print(f"Error processing reaction SMILES at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_cn_formation
