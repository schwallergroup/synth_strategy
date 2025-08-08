#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects if the synthesis involves multiple N-alkylation or N-acylation steps.
    """
    n_functionalization_count = 0

    def dfs_traverse(node):
        nonlocal n_functionalization_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-alkylation (NH to N-C) or N-acylation (NH to N-C=O)
            nh_pattern = Chem.MolFromSmarts("[#7H]")
            n_alkyl_pattern = Chem.MolFromSmarts("[#7]-[#6]")
            n_acyl_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and any(
                r and r.HasSubstructMatch(nh_pattern) for r in reactant_mols if r
            ):
                if product_mol.HasSubstructMatch(n_alkyl_pattern) or product_mol.HasSubstructMatch(
                    n_acyl_pattern
                ):
                    print(
                        f"Detected N-functionalization step at depth {node.get('depth', 'unknown')}"
                    )
                    n_functionalization_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_functionalization_count >= 2  # Return True if at least 2 N-functionalization steps
