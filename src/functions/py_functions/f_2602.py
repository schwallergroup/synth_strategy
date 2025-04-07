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
    This function detects a strategy involving sequential functionalization
    of an aromatic ring with electron-withdrawing and electron-donating groups.
    """
    # Track functional group installations
    functional_group_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for nitro group installation
                nitro_pattern = Chem.MolFromSmarts("[#6]~[#7](~[#8])(~[#8])")
                product_nitro_count = len(product.GetSubstructMatches(nitro_pattern))
                reactants_nitro_count = sum(
                    len(r.GetSubstructMatches(nitro_pattern)) for r in reactants
                )

                if product_nitro_count > reactants_nitro_count:
                    functional_group_sequence.append(("nitro", depth))
                    print(f"Detected nitro group installation at depth {depth}")

                # Check for methoxy group installation
                methoxy_pattern = Chem.MolFromSmarts("[#6]~[#8]~[#6&D1]")
                product_methoxy_count = len(
                    product.GetSubstructMatches(methoxy_pattern)
                )
                reactants_methoxy_count = sum(
                    len(r.GetSubstructMatches(methoxy_pattern)) for r in reactants
                )

                if product_methoxy_count > reactants_methoxy_count:
                    functional_group_sequence.append(("methoxy", depth))
                    print(f"Detected methoxy group installation at depth {depth}")

                # Check for esterification
                acid_pattern = Chem.MolFromSmarts("[#6](=[#8])([#8;H1])")
                ester_pattern = Chem.MolFromSmarts("[#6](=[#8])([#8][#6&D1])")

                product_acid_count = len(product.GetSubstructMatches(acid_pattern))
                reactants_acid_count = sum(
                    len(r.GetSubstructMatches(acid_pattern)) for r in reactants
                )
                product_ester_count = len(product.GetSubstructMatches(ester_pattern))
                reactants_ester_count = sum(
                    len(r.GetSubstructMatches(ester_pattern)) for r in reactants
                )

                if (
                    product_ester_count > reactants_ester_count
                    and reactants_acid_count > product_acid_count
                ):
                    functional_group_sequence.append(("ester", depth))
                    print(f"Detected esterification at depth {depth}")

                # Check for benzylic functionalization
                benzylic_pattern = Chem.MolFromSmarts("[#6;a]-[#6;!a]")
                benzylic_br_pattern = Chem.MolFromSmarts("[#6;a]-[#6;!a]-[#35]")

                if (
                    any(
                        len(r.GetSubstructMatches(benzylic_br_pattern)) > 0
                        for r in reactants
                    )
                    and len(product.GetSubstructMatches(benzylic_pattern)) > 0
                ):
                    functional_group_sequence.append(("benzylic", depth))
                    print(f"Detected benzylic functionalization at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth (higher depth = earlier in synthesis)
    functional_group_sequence.sort(key=lambda x: -x[1])

    # Check if we have at least 3 sequential functionalizations
    strategy_present = len(functional_group_sequence) >= 3

    # Check if electron-withdrawing groups (nitro) are installed before electron-donating groups (methoxy)
    if strategy_present:
        group_types = [g[0] for g in functional_group_sequence]
        nitro_index = group_types.index("nitro") if "nitro" in group_types else -1
        methoxy_index = group_types.index("methoxy") if "methoxy" in group_types else -1

        if nitro_index != -1 and methoxy_index != -1:
            strategy_present = (
                nitro_index < methoxy_index
            )  # Lower index means earlier in the sequence

    print(f"Sequential functionalization strategy detected: {strategy_present}")
    print(f"Functionalization sequence: {functional_group_sequence}")
    return strategy_present
