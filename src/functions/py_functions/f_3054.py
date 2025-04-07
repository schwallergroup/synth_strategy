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
    This function detects a synthetic strategy that employs benzyl groups for protection/functionalization
    of both phenolic OH and carboxylic acid groups, with at least two benzyl protection events.
    """
    # Track benzyl protection events
    benzyl_protection_events = []

    # SMARTS patterns for detecting relevant structures
    benzyl_ester_patt = Chem.MolFromSmarts("[#6]C(=O)OC[c]1[cH][cH][cH][cH][cH]1")
    benzyl_ether_patt = Chem.MolFromSmarts("[#8]C[c]1[cH][cH][cH][cH][cH]1")

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles)

                # Check if product contains benzyl ester or benzyl ether
                has_benzyl_ester = (
                    product.HasSubstructMatch(benzyl_ester_patt) if product else False
                )
                has_benzyl_ether = (
                    product.HasSubstructMatch(benzyl_ether_patt) if product else False
                )

                # Check if any reactant has phenol or carboxylic acid
                has_phenol = any(
                    r and r.HasSubstructMatch(Chem.MolFromSmarts("[OH][c]"))
                    for r in reactants
                )
                has_carboxylic_acid = any(
                    r and r.HasSubstructMatch(Chem.MolFromSmarts("[#6]C(=O)[OH]"))
                    for r in reactants
                )

                # Detect benzyl protection events
                if (has_benzyl_ester and has_carboxylic_acid) or (
                    has_benzyl_ether and has_phenol
                ):
                    benzyl_protection_events.append(
                        {
                            "depth": depth,
                            "type": "ester" if has_benzyl_ester else "ether",
                        }
                    )
                    print(
                        f"Detected benzyl {'ester' if has_benzyl_ester else 'ether'} formation at depth {depth}"
                    )
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least two benzyl protection events
    if len(benzyl_protection_events) >= 2:
        # Check if we have both types (ester and ether)
        has_ester = any(event["type"] == "ester" for event in benzyl_protection_events)
        has_ether = any(event["type"] == "ether" for event in benzyl_protection_events)

        # Check if protections occur at different stages (depths)
        depths = set(event["depth"] for event in benzyl_protection_events)

        print(
            f"Found {len(benzyl_protection_events)} benzyl protection events at depths {depths}"
        )

        # Strategy is present if we have multiple benzyl protections at different depths
        return len(depths) >= 2

    return False
