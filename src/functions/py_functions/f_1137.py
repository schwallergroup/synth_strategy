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
    Detects a synthetic strategy where a trifluoromethoxy group is preserved throughout
    the synthesis and culminates in isocyanate formation
    """
    # Initialize tracking variables
    has_trifluoromethoxy_in_final = False
    has_isocyanate_formation = False

    # SMARTS patterns
    trifluoromethoxy_pattern = Chem.MolFromSmarts("[#8]-[#6](-[#9])(-[#9])-[#9]")
    amine_pattern = Chem.MolFromSmarts("[NH2]-c")
    isocyanate_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#8]")

    # Track the final product
    final_product = None

    def dfs_traverse(node, is_root=True):
        nonlocal has_trifluoromethoxy_in_final, has_isocyanate_formation, final_product

        if is_root and node["type"] == "mol":
            # This is the final product
            final_product = Chem.MolFromSmiles(node["smiles"])
            if final_product:
                if final_product.HasSubstructMatch(trifluoromethoxy_pattern):
                    has_trifluoromethoxy_in_final = True
                    print("Detected trifluoromethoxy group in final product")

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants = parts[0].split(".")
                    product = parts[2]

                    # Convert to RDKit molecules
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    # Check for isocyanate formation
                    if product_mol and product_mol.HasSubstructMatch(
                        isocyanate_pattern
                    ):
                        if any(
                            r and r.HasSubstructMatch(amine_pattern)
                            for r in reactant_mols
                            if r
                        ):
                            has_isocyanate_formation = True
                            print("Detected isocyanate formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, is_root=False)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = has_trifluoromethoxy_in_final and has_isocyanate_formation

    if strategy_present:
        print("Detected preserved trifluoromethoxy with isocyanate formation strategy")

    return strategy_present
