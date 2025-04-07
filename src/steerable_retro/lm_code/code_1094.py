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
    This function detects if the synthetic route involves transformations
    of aromatic C-X bonds to different C-Y bonds (e.g., C-F to C-O, C-Br to C-B).
    """
    c_x_to_c_y_found = False

    def dfs_traverse(node):
        nonlocal c_x_to_c_y_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                # Check for aromatic C-X bonds in reactants
                c_f_pattern = Chem.MolFromSmarts("[c]-[#9]")
                c_cl_pattern = Chem.MolFromSmarts("[c]-[#17]")
                c_br_pattern = Chem.MolFromSmarts("[c]-[#35]")
                c_i_pattern = Chem.MolFromSmarts("[c]-[#53]")

                # Check for aromatic C-Y bonds in product
                c_o_pattern = Chem.MolFromSmarts("[c]-[#8]")
                c_n_pattern = Chem.MolFromSmarts("[c]-[#7]")
                c_b_pattern = Chem.MolFromSmarts("[c]-[#5]")
                c_c_pattern = Chem.MolFromSmarts("[c]-[#6]")

                reactants_have_c_x = any(
                    r
                    and (
                        r.HasSubstructMatch(c_f_pattern)
                        or r.HasSubstructMatch(c_cl_pattern)
                        or r.HasSubstructMatch(c_br_pattern)
                        or r.HasSubstructMatch(c_i_pattern)
                    )
                    for r in reactant_mols
                    if r
                )

                if product_mol and reactants_have_c_x:
                    product_has_c_y = (
                        product_mol.HasSubstructMatch(c_o_pattern)
                        or product_mol.HasSubstructMatch(c_n_pattern)
                        or product_mol.HasSubstructMatch(c_b_pattern)
                        or product_mol.HasSubstructMatch(c_c_pattern)
                    )

                    if product_has_c_y:
                        print(f"Aromatic C-X to C-Y transformation detected in reaction: {rsmi}")
                        c_x_to_c_y_found = True
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return c_x_to_c_y_found
