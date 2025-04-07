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
    This function detects a synthetic strategy involving multiple protection/deprotection steps.
    """
    protection_count = 0
    deprotection_count = 0

    # SMARTS patterns for common protecting groups
    boc_pattern = Chem.MolFromSmarts("[NX3]C(=O)OC(C)(C)C")
    silyl_pattern = Chem.MolFromSmarts("[OX2][Si]")
    nosyl_pattern = Chem.MolFromSmarts("[NX3][S](=O)(=O)c1ccc([N+](=O)[O-])cc1")

    def dfs_traverse(node):
        nonlocal protection_count, deprotection_count

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactants_mol and product_mol:
                        # Check for protection (appearance of protecting group)
                        if (
                            (
                                not reactants_mol.HasSubstructMatch(boc_pattern)
                                and product_mol.HasSubstructMatch(boc_pattern)
                            )
                            or (
                                not reactants_mol.HasSubstructMatch(silyl_pattern)
                                and product_mol.HasSubstructMatch(silyl_pattern)
                            )
                            or (
                                not reactants_mol.HasSubstructMatch(nosyl_pattern)
                                and product_mol.HasSubstructMatch(nosyl_pattern)
                            )
                        ):
                            protection_count += 1
                            print(f"Protection step detected: {rsmi}")

                        # Check for deprotection (disappearance of protecting group)
                        if (
                            (
                                reactants_mol.HasSubstructMatch(boc_pattern)
                                and not product_mol.HasSubstructMatch(boc_pattern)
                            )
                            or (
                                reactants_mol.HasSubstructMatch(silyl_pattern)
                                and not product_mol.HasSubstructMatch(silyl_pattern)
                            )
                            or (
                                reactants_mol.HasSubstructMatch(nosyl_pattern)
                                and not product_mol.HasSubstructMatch(nosyl_pattern)
                            )
                        ):
                            deprotection_count += 1
                            print(f"Deprotection step detected: {rsmi}")
                except:
                    print("Error processing SMILES in multiple_protection_deprotection_strategy")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if there are at least 2 protection/deprotection steps combined
    return (protection_count + deprotection_count) >= 2
