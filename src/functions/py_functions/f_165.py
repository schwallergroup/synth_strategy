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
    This function detects a synthetic strategy involving chain extension through
    etherification (C-O-C bond formation).
    """
    etherification_found = False

    def dfs_traverse(node):
        nonlocal etherification_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                try:
                    reactants = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")
                    ]
                    product = Chem.MolFromSmiles(products_smiles)

                    # Check for etherification: one reactant has OH, another has leaving group
                    if product and len(reactants) >= 2:
                        oh_pattern = Chem.MolFromSmarts("[OH]")
                        br_pattern = Chem.MolFromSmarts("[Br]")

                        has_oh = any(
                            r and r.HasSubstructMatch(oh_pattern)
                            for r in reactants
                            if r
                        )
                        has_br = any(
                            r and r.HasSubstructMatch(br_pattern)
                            for r in reactants
                            if r
                        )

                        # Check if product has new ether linkage
                        ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
                        if (
                            has_oh
                            and has_br
                            and product.HasSubstructMatch(ether_pattern)
                        ):
                            # Count ethers in reactants vs product
                            reactant_ether_count = sum(
                                len(r.GetSubstructMatches(ether_pattern))
                                for r in reactants
                                if r and r.HasSubstructMatch(ether_pattern)
                            )
                            product_ether_count = len(
                                product.GetSubstructMatches(ether_pattern)
                            )

                            if product_ether_count > reactant_ether_count:
                                etherification_found = True
                                print("Etherification for chain extension detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return etherification_found
