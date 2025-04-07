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
    This function detects a synthetic strategy involving multiple oxygen-containing
    functional group manipulations (alcohols, ethers, esters, carboxylic acids).
    """
    o_fg_transformations = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal o_fg_transformations, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
            ]
            product = Chem.MolFromSmiles(product_smiles)

            if not product or not reactants:
                return

            # Check for ester reduction
            ester_patt = Chem.MolFromSmarts("[C](=[O])[O][C]")
            alcohol_patt = Chem.MolFromSmarts("[O]-[CH2]")

            if any(r.HasSubstructMatch(ester_patt) for r in reactants):
                if product.HasSubstructMatch(alcohol_patt):
                    o_fg_transformations += 1
                    print(f"Ester reduction detected in reaction: {rsmi}")

            # Check for ester hydrolysis
            acid_patt = Chem.MolFromSmarts("[C](=[O])[O][H]")

            if any(r.HasSubstructMatch(ester_patt) for r in reactants):
                if product.HasSubstructMatch(acid_patt):
                    o_fg_transformations += 1
                    print(f"Ester hydrolysis detected in reaction: {rsmi}")

            # Check for O-alkylation
            phenol_patt = Chem.MolFromSmarts("c-[O][H]")
            benzyl_ether_patt = Chem.MolFromSmarts("c-[O]-[CH2]-c")

            if any(r.HasSubstructMatch(phenol_patt) for r in reactants):
                if product.HasSubstructMatch(benzyl_ether_patt):
                    o_fg_transformations += 1
                    print(f"O-alkylation detected in reaction: {rsmi}")

            # Check for ether cleavage
            ether_patt = Chem.MolFromSmarts("[O]-[CH3]")

            if any(r.HasSubstructMatch(ether_patt) for r in reactants):
                if product.HasSubstructMatch(phenol_patt):
                    o_fg_transformations += 1
                    print(f"Ether cleavage detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Strategy criteria: at least 3 oxygen functional group transformations
    has_strategy = o_fg_transformations >= 3
    print(f"Oxygen functional group transformations: {o_fg_transformations}")
    print(f"Strategy detected: {has_strategy}")

    return has_strategy
