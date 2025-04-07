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
    Detects a synthetic strategy involving sequential transformations:
    thioether -> sulfone -> nitro -> amine -> isocyanate
    with preserved trifluoromethoxy group
    """
    # Initialize tracking variables
    has_trifluoromethoxy = False
    has_thioether_to_sulfone = False
    has_nitro_to_amine = False
    has_amine_to_isocyanate = False

    # SMARTS patterns
    trifluoromethoxy_pattern = Chem.MolFromSmarts("[#8]-[#6](-[#9])(-[#9])-[#9]")
    thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-c")
    sulfone_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])(=[#8])-c")
    nitro_pattern = Chem.MolFromSmarts("[$([#7](=[#8])~[#8])]")
    amine_pattern = Chem.MolFromSmarts("[NH2]-c")
    isocyanate_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#8]")

    def dfs_traverse(node):
        nonlocal has_trifluoromethoxy, has_thioether_to_sulfone, has_nitro_to_amine, has_amine_to_isocyanate

        if node["type"] == "mol":
            # Check for trifluoromethoxy group in molecules
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(trifluoromethoxy_pattern):
                has_trifluoromethoxy = True

        elif node["type"] == "reaction":
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

                    # Check for thioether to sulfone transformation
                    if product_mol and product_mol.HasSubstructMatch(sulfone_pattern):
                        if any(
                            r and r.HasSubstructMatch(thioether_pattern) for r in reactant_mols if r
                        ):
                            has_thioether_to_sulfone = True
                            print("Detected thioether to sulfone transformation")

                    # Check for nitro to amine transformation
                    if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                        if any(
                            r and r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r
                        ):
                            has_nitro_to_amine = True
                            print("Detected nitro to amine transformation")

                    # Check for amine to isocyanate transformation
                    if product_mol and product_mol.HasSubstructMatch(isocyanate_pattern):
                        if any(
                            r and r.HasSubstructMatch(amine_pattern) for r in reactant_mols if r
                        ):
                            has_amine_to_isocyanate = True
                            print("Detected amine to isocyanate transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if all required transformations are present
    strategy_present = (
        has_trifluoromethoxy
        and has_thioether_to_sulfone
        and has_nitro_to_amine
        and has_amine_to_isocyanate
    )

    if strategy_present:
        print(
            "Detected complete sulfone-nitro-amine-isocyanate sequence with preserved trifluoromethoxy group"
        )

    return strategy_present
