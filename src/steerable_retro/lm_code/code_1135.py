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
    Detects a strategy involving early S-methylation followed by
    later stage aromatic functionalization (nitration, reduction, isocyanate formation)
    """
    # Initialize tracking variables
    has_s_methylation = False
    s_methylation_depth = -1
    has_nitration = False
    nitration_depth = -1
    has_amine_formation = False
    amine_formation_depth = -1
    has_isocyanate_formation = False
    isocyanate_formation_depth = -1

    # SMARTS patterns
    thiol_pattern = Chem.MolFromSmarts("[#16H]-c")
    methylthio_pattern = Chem.MolFromSmarts("[#6]-[#16]-c")
    nitro_pattern = Chem.MolFromSmarts("[$([#7](=[#8])~[#8])]")
    amine_pattern = Chem.MolFromSmarts("[NH2]-c")
    isocyanate_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#8]")

    def dfs_traverse(node, depth=0):
        nonlocal has_s_methylation, s_methylation_depth, has_nitration, nitration_depth
        nonlocal has_amine_formation, amine_formation_depth, has_isocyanate_formation, isocyanate_formation_depth

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

                    # Check for S-methylation
                    if product_mol and product_mol.HasSubstructMatch(methylthio_pattern):
                        if any(
                            r and r.HasSubstructMatch(thiol_pattern) for r in reactant_mols if r
                        ):
                            has_s_methylation = True
                            s_methylation_depth = depth
                            print(f"Detected S-methylation at depth {depth}")

                    # Check for nitration
                    if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                        has_nitration = True
                        nitration_depth = depth
                        print(f"Detected nitration at depth {depth}")

                    # Check for amine formation
                    if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                        if any(
                            r and r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r
                        ):
                            has_amine_formation = True
                            amine_formation_depth = depth
                            print(f"Detected amine formation at depth {depth}")

                    # Check for isocyanate formation
                    if product_mol and product_mol.HasSubstructMatch(isocyanate_pattern):
                        if any(
                            r and r.HasSubstructMatch(amine_pattern) for r in reactant_mols if r
                        ):
                            has_isocyanate_formation = True
                            isocyanate_formation_depth = depth
                            print(f"Detected isocyanate formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present (early S-methylation followed by late functionalization)
    strategy_present = (
        has_s_methylation
        and has_nitration
        and has_amine_formation
        and has_isocyanate_formation
        and s_methylation_depth
        > nitration_depth
        > amine_formation_depth
        > isocyanate_formation_depth
    )

    if strategy_present:
        print("Detected early S-methylation with late-stage functionalization strategy")

    return strategy_present
