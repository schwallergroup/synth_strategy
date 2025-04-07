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
    This function detects if the synthesis route involves a protection-deprotection sequence,
    specifically looking for carboxylic acid → tert-butyl ester → carboxylic acid pattern.
    """
    protection_steps = []
    deprotection_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if all(reactant_mols) and product_mol:
                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#8]-[#6](=[#8])-[#6]")
                # Check for tert-butyl ester in product
                tbutyl_ester_pattern = Chem.MolFromSmarts(
                    "[#6]C([#6])([#6])[#6]-[#8]-[#6](=[#8])-[#6]"
                )

                # Protection: acid → tert-butyl ester
                if any(
                    mol.HasSubstructMatch(carboxylic_acid_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(tbutyl_ester_pattern):
                    protection_steps.append(depth)
                    print(f"Detected protection (acid → tert-butyl ester) at depth {depth}")

                # Deprotection: tert-butyl ester → acid
                if any(
                    mol.HasSubstructMatch(tbutyl_ester_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                    deprotection_steps.append(depth)
                    print(f"Detected deprotection (tert-butyl ester → acid) at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if both protection and deprotection occurred
    has_protection_deprotection = len(protection_steps) > 0 and len(deprotection_steps) > 0

    return has_protection_deprotection
