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
    Detects if the synthesis route uses a late-stage SNAr coupling (depth 0-1)
    between an amine nucleophile and a halogen-substituted aromatic/heteroaromatic.
    """
    snar_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_detected

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine nucleophile
                amine_pattern = Chem.MolFromSmarts("[NH1,NH2]")
                # Check for halogen-substituted aromatic/heteroaromatic
                halo_aromatic_pattern = Chem.MolFromSmarts("c[Cl,F,Br]")
                halo_hetero_pattern = Chem.MolFromSmarts("n[Cl,F,Br]")

                # Check reactants for patterns
                amine_present = False
                halo_present = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(amine_pattern):
                            amine_present = True
                        if mol.HasSubstructMatch(halo_aromatic_pattern) or mol.HasSubstructMatch(
                            halo_hetero_pattern
                        ):
                            halo_present = True

                # Check if product has new C-N bond
                if amine_present and halo_present:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        # If both patterns were in reactants and product has fewer halogens,
                        # likely an SNAr reaction occurred
                        reactant_halo_count = sum(
                            len(
                                Chem.MolFromSmiles(r).GetSubstructMatches(
                                    Chem.MolFromSmarts("[Cl,F,Br]")
                                )
                            )
                            for r in reactants
                            if Chem.MolFromSmiles(r)
                        )
                        product_halo_count = len(
                            prod_mol.GetSubstructMatches(Chem.MolFromSmarts("[Cl,F,Br]"))
                        )

                        if product_halo_count < reactant_halo_count:
                            print(f"SNAr coupling detected at depth {depth}")
                            snar_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return snar_detected
