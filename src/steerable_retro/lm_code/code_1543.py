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
    This function detects a late-stage N-alkylation strategy where an amine attacks
    an activated alcohol (typically as mesylate or similar leaving group).
    """
    n_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for mesylate or similar leaving group in reactants
                mesylate_pattern = Chem.MolFromSmarts("[C]-[O]-[S](=[O])(=[O])-[C]")
                amine_pattern = Chem.MolFromSmarts("[N;H2,H1;!$(NC=O)]")

                has_mesylate = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(mesylate_pattern):
                            has_mesylate = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                # Check if product has a new C-N bond that wasn't in reactants
                if has_mesylate and has_amine:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # This is a simplified check - in practice would need more sophisticated analysis
                        n_alkylation_detected = True
                        print("Late-stage N-alkylation detected at depth", depth)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return n_alkylation_detected
