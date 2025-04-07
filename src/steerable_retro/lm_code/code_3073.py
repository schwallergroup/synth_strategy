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
    Detects a sequence of nitrogen functionalization steps:
    primary amine → secondary amine → tertiary amine/amide
    """
    # Track nitrogen functionalization steps
    primary_to_secondary = False
    secondary_to_tertiary = False

    def dfs_traverse(node):
        nonlocal primary_to_secondary, secondary_to_tertiary

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Try to create molecules
            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if None in reactant_mols or product_mol is None:
                    return

                # Check for primary amine → secondary amine transformation
                primary_amine_pattern = Chem.MolFromSmarts("[#7;H2]")
                secondary_amine_pattern = Chem.MolFromSmarts("[#7;H1;!$(N-C=O)]")

                has_primary_amine = any(
                    mol.HasSubstructMatch(primary_amine_pattern) for mol in reactant_mols
                )
                has_secondary_amine_product = product_mol.HasSubstructMatch(secondary_amine_pattern)

                if has_primary_amine and has_secondary_amine_product:
                    primary_to_secondary = True
                    print("Found primary amine → secondary amine transformation")

                # Check for secondary amine → tertiary amine/amide transformation
                tertiary_nitrogen_pattern = Chem.MolFromSmarts("[#7;H0]")

                has_secondary_amine = any(
                    mol.HasSubstructMatch(secondary_amine_pattern) for mol in reactant_mols
                )
                has_tertiary_nitrogen = product_mol.HasSubstructMatch(tertiary_nitrogen_pattern)

                if has_secondary_amine and has_tertiary_nitrogen:
                    secondary_to_tertiary = True
                    print("Found secondary amine → tertiary nitrogen transformation")

            except:
                pass  # Skip if there's an error processing the molecules

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the complete sequence
    result = primary_to_secondary and secondary_to_tertiary

    print(f"Nitrogen functionalization sequence detected: {result}")
    return result
