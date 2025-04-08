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
    This function detects peptide synthesis strategy using NHS activation.
    It looks for:
    1. Multiple amide bond formations
    2. Use of NHS activation in at least one step
    3. Maintenance of an allyl ether protecting group
    """
    amide_formations = 0
    nhs_activation_used = False
    allyl_ether_maintained = False

    # SMARTS patterns
    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
    nhs_pattern = Chem.MolFromSmarts("[O][N]1[C](=[O])[CH2][CH2][C]1=[O]")
    allyl_ether_pattern = Chem.MolFromSmarts("[CH2]=[CH][CH2][O]")

    def dfs_traverse(node):
        nonlocal amide_formations, nhs_activation_used, allyl_ether_maintained

        if node["type"] == "mol":
            # Check if final product has allyl ether group
            if node.get("in_stock", False) == False:  # Not a starting material
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(allyl_ether_pattern):
                    allyl_ether_maintained = True
                    print(f"Allyl ether found in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for NHS in reactants
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(nhs_pattern):
                        nhs_activation_used = True
                        print(f"NHS activation detected in reaction: {rsmi}")

                # Count amide formations by comparing reactants and products
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    product_amide_count = len(product_mol.GetSubstructMatches(amide_pattern))

                    reactant_amide_count = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            reactant_amide_count += len(mol.GetSubstructMatches(amide_pattern))

                    # If product has more amides than reactants combined, amide formation occurred
                    if product_amide_count > reactant_amide_count:
                        amide_formations += 1
                        print(f"Amide formation detected in reaction: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Criteria for peptide synthesis with NHS activation
    result = amide_formations >= 2 and nhs_activation_used and allyl_ether_maintained
    print(f"Peptide synthesis with NHS activation: {result}")
    print(
        f"Amide formations: {amide_formations}, NHS used: {nhs_activation_used}, Allyl ether maintained: {allyl_ether_maintained}"
    )

    return result
