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
    This function detects a synthetic strategy involving protection/deprotection sequences.
    """
    protection_count = 0

    def dfs_traverse(node):
        nonlocal protection_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for common protection reactions

                    # 1. Boc protection
                    boc_pattern = Chem.MolFromSmarts(
                        "[#6]-[#6](-[#6])(-[#6])-[#8]-[#6](=[#8])-[#7]"
                    )

                    # 2. Ester formation (carboxylic acid protection)
                    ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])")

                    # 3. Benzyl ether formation (alcohol/phenol protection)
                    benzyl_pattern = Chem.MolFromSmarts("c1ccccc1-[#6]-[#8]")

                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    if product_mol:
                        # Check if product has protection group but not all reactants do
                        for pattern, name in [
                            (boc_pattern, "Boc"),
                            (ester_pattern, "Ester"),
                            (benzyl_pattern, "Benzyl"),
                        ]:
                            if product_mol.HasSubstructMatch(pattern):
                                all_reactants_have_pattern = all(
                                    r and r.HasSubstructMatch(pattern)
                                    for r in reactant_mols
                                )
                                if not all_reactants_have_pattern:
                                    protection_count += 1
                                    print(
                                        f"{name} protection detected in reaction: {rsmi}"
                                    )
                except:
                    pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total protection reactions detected: {protection_count}")
    return protection_count >= 2
