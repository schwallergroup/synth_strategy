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
    Detects a strategy involving alcohol activation (via mesylation) followed by N-alkylation,
    with preceding functional group transformations from carboxylic acid to alcohol via ester.
    """
    # Track if we've found each step in the strategy
    found_esterification = False
    found_ester_reduction = False
    found_mesylation = False
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_esterification, found_ester_reduction, found_mesylation, found_n_alkylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and all(r for r in reactants):
                    # Check for esterification (carboxylic acid + alcohol → ester)
                    acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8;H1]")
                    ester_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8]-[#6]")
                    methanol_pattern = Chem.MolFromSmarts("[#8;H1]-[#6;H3]")

                    has_acid = any(r.HasSubstructMatch(acid_pattern) for r in reactants)
                    has_methanol = any(r.HasSubstructMatch(methanol_pattern) for r in reactants)
                    has_ester_product = product.HasSubstructMatch(ester_pattern)

                    if has_acid and has_methanol and has_ester_product:
                        print("Found esterification step")
                        found_esterification = True

                    # Check for ester reduction (ester → alcohol)
                    if any(r.HasSubstructMatch(ester_pattern) for r in reactants):
                        alcohol_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#8;H1]")
                        if product.HasSubstructMatch(alcohol_pattern):
                            print("Found ester reduction step")
                            found_ester_reduction = True

                    # Check for mesylation (alcohol → mesylate)
                    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                    mesylate_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#16](=[#8])(=[#8])-[#6]")
                    mesyl_chloride_pattern = Chem.MolFromSmarts("[Cl]-[#16](=[#8])(=[#8])-[#6]")

                    has_alcohol = any(r.HasSubstructMatch(alcohol_pattern) for r in reactants)
                    has_mesyl_chloride = any(
                        r.HasSubstructMatch(mesyl_chloride_pattern) for r in reactants
                    )
                    has_mesylate_product = product.HasSubstructMatch(mesylate_pattern)

                    if has_alcohol and has_mesyl_chloride and has_mesylate_product:
                        print("Found mesylation step")
                        found_mesylation = True

                    # Check for N-alkylation (secondary amine + mesylate → tertiary amine)
                    sec_amine_pattern = Chem.MolFromSmarts("[#6]-[#7;H1]-[#6]")
                    tert_amine_pattern = Chem.MolFromSmarts("[#6]-[#7](-[#6])-[#6]")

                    has_sec_amine = any(r.HasSubstructMatch(sec_amine_pattern) for r in reactants)
                    has_mesylate = any(r.HasSubstructMatch(mesylate_pattern) for r in reactants)
                    has_tert_amine_product = product.HasSubstructMatch(tert_amine_pattern)

                    if has_sec_amine and has_mesylate and has_tert_amine_product:
                        print("Found N-alkylation step")
                        found_n_alkylation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # The strategy is present if we found all the key steps
    strategy_present = (
        found_esterification and found_ester_reduction and found_mesylation and found_n_alkylation
    )
    print(f"Alcohol activation N-alkylation strategy detected: {strategy_present}")
    return strategy_present
