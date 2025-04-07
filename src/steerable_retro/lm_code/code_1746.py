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
    Detects late-stage SNAr coupling (final step) of two complex fragments.
    """
    snar_coupling_found = False

    def dfs_traverse(node):
        nonlocal snar_coupling_found

        # We're looking for the final step (depth 0)
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a final step (no parent reaction)
            if product:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and len(reactant_mols) >= 2:
                    # Look for halogen on heteroaromatic in reactants
                    halo_heteroaromatic = Chem.MolFromSmarts("[n]([c][Cl,Br,I,F])")
                    # Look for amine nucleophile in reactants
                    amine_pattern = Chem.MolFromSmarts("[NH2][c]")

                    has_halo = any(
                        mol.HasSubstructMatch(halo_heteroaromatic) for mol in reactant_mols if mol
                    )
                    has_amine = any(
                        mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                    )

                    # Look for C-N bond formation in product
                    cn_bond_pattern = Chem.MolFromSmarts("[n]([c][N][c])")
                    has_cn_bond = (
                        product_mol.HasSubstructMatch(cn_bond_pattern) if product_mol else False
                    )

                    if has_halo and has_amine and has_cn_bond:
                        # Check if this is a complex coupling (both fragments have significant complexity)
                        complex_fragments = sum(
                            1 for mol in reactant_mols if mol and mol.GetNumAtoms() > 8
                        )

                        if complex_fragments >= 2:
                            snar_coupling_found = True
                            print("Found late-stage SNAr coupling of complex fragments")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return snar_coupling_found
