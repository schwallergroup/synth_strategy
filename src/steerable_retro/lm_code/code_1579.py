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
    Detects if the synthesis involves a late-stage N-benzylation of a piperazine or similar nitrogen heterocycle.
    Late stage means it occurs in the first half of the synthesis (low depth in retrosynthetic tree).
    """
    benzylation_detected = False
    max_depth = 0
    benzylation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal benzylation_detected, max_depth, benzylation_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an N-benzylation reaction
            try:
                prod_mol = Chem.MolFromSmiles(product)
                react_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Look for benzyl group in product
                benzyl_patt = Chem.MolFromSmarts("[c][CH2][N;!$(NC=O)]")
                if prod_mol and prod_mol.HasSubstructMatch(benzyl_patt):
                    # Check if benzyl group was added in this step
                    benzyl_in_reactants = any(
                        m and m.HasSubstructMatch(benzyl_patt) for m in react_mols if m
                    )

                    if not benzyl_in_reactants:
                        print(f"N-benzylation detected at depth {depth}")
                        benzylation_detected = True
                        benzylation_depth = depth
            except:
                print("Error in processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if benzylation is in the first half of the synthesis (late stage)
    if benzylation_detected and benzylation_depth is not None:
        is_late_stage = benzylation_depth <= max_depth / 2
        print(
            f"Benzylation depth: {benzylation_depth}, Max depth: {max_depth}, Is late stage: {is_late_stage}"
        )
        return is_late_stage

    return False
