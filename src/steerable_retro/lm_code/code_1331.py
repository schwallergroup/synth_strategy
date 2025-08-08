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
    This function detects a synthetic strategy involving both oxidation and reduction steps
    in a specific sequence (e.g., alcohol oxidation followed by oxime reduction).
    """
    # Track oxidation and reduction reactions
    oxidation_reactions = []
    reduction_reactions = []
    oxidation_depths = []
    reduction_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and any(r for r in reactant_mols):
                    # Check for alcohol to aldehyde/ketone/acid oxidation
                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[CH2]-[OH]"))
                        for r in reactant_mols
                        if r
                    ) and (
                        product_mol.HasSubstructMatch(Chem.MolFromSmarts("[CH]=O"))
                        or product_mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)-O"))
                    ):
                        oxidation_reactions.append(rsmi)
                        oxidation_depths.append(depth)
                        print(f"Found oxidation reaction at depth {depth}: {rsmi}")

                    # Check for oxime/imine reduction
                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[CH]=[N]"))
                        for r in reactant_mols
                        if r
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]-[N]")):
                        reduction_reactions.append(rsmi)
                        reduction_depths.append(depth)
                        print(f"Found reduction reaction at depth {depth}: {rsmi}")

                    # Check for aldehyde/ketone reduction
                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[CH]=O"))
                        for r in reactant_mols
                        if r
                    ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]-[OH]")):
                        reduction_reactions.append(rsmi)
                        reduction_depths.append(depth)
                        print(f"Found reduction reaction at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both oxidation and reduction reactions
    has_oxidation = len(oxidation_reactions) > 0
    has_reduction = len(reduction_reactions) > 0

    # Check if oxidation occurs before reduction in the synthetic direction
    # (which means higher depth in retrosynthetic direction)
    oxidation_before_reduction = False
    if has_oxidation and has_reduction:
        max_oxidation_depth = max(oxidation_depths) if oxidation_depths else -1
        min_reduction_depth = min(reduction_depths) if reduction_depths else float("inf")
        oxidation_before_reduction = max_oxidation_depth > min_reduction_depth

    strategy_present = has_oxidation and has_reduction and oxidation_before_reduction
    print(f"Oxidation-reduction sequence strategy detected: {strategy_present}")
    print(
        f"Found {len(oxidation_reactions)} oxidation reactions and {len(reduction_reactions)} reduction reactions"
    )

    return strategy_present
