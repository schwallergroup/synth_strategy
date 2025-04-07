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
    Detects if the synthesis route involves multiple amide coupling reactions.
    """
    amide_coupling_count = 0

    def dfs_traverse(node):
        nonlocal amide_coupling_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide coupling: product has C(=O)N pattern and reactants have C(=O)O and NH2/NH patterns
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[N]")):
                # Check reactants for carboxylic acid and amine
                has_carboxylic_acid = False
                has_amine = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O]")):
                            has_carboxylic_acid = True
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[N]")):
                            has_amine = True

                if has_carboxylic_acid and has_amine:
                    amide_coupling_count += 1
                    print(
                        f"Detected amide coupling at depth {node['metadata'].get('depth', 'unknown')}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return amide_coupling_count >= 2  # Return True if at least 2 amide couplings are found
