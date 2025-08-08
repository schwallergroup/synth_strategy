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
    This function detects a linear synthesis strategy using malononitrile as a key building block,
    with a sequence of Knoevenagel condensation, reduction, and alkylation.
    """
    # Initialize flags for each step in the strategy
    has_knoevenagel = False
    has_reduction_after_knoevenagel = False
    has_alkylation_after_reduction = False

    # Track the depth of reactions for sequence analysis
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal has_knoevenagel, has_reduction_after_knoevenagel, has_alkylation_after_reduction

        if node["type"] == "reaction":
            # Store reaction information with depth
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
            product = Chem.MolFromSmiles(product_str) if product_str else None

            if product and reactants:
                reaction_sequence.append((depth, reactants, product, rsmi))

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Traverse the route
    dfs_traverse(route)

    # Sort reactions by depth (higher depth = earlier in synthesis)
    reaction_sequence.sort(key=lambda x: x[0], reverse=True)

    # Check for the specific sequence pattern
    for i, (depth, reactants, product, rsmi) in enumerate(reaction_sequence):
        # Check for Knoevenagel condensation (aldehyde + malononitrile → α,β-unsaturated nitrile)
        if not has_knoevenagel:
            aldehyde_pattern = Chem.MolFromSmarts("[#6][CH]=O")
            malononitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])")
            unsaturated_nitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])=[#6]")

            has_aldehyde = any(mol.HasSubstructMatch(aldehyde_pattern) for mol in reactants if mol)
            has_malononitrile = any(
                mol.HasSubstructMatch(malononitrile_pattern) for mol in reactants if mol
            )
            product_has_unsaturated_nitrile = product and product.HasSubstructMatch(
                unsaturated_nitrile_pattern
            )

            if has_aldehyde and has_malononitrile and product_has_unsaturated_nitrile:
                has_knoevenagel = True
                print(f"Found Knoevenagel condensation at depth {depth}")
                continue

        # Check for reduction of α,β-unsaturated nitrile after Knoevenagel
        if has_knoevenagel and not has_reduction_after_knoevenagel and i > 0:
            unsaturated_nitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])=[#6]")
            saturated_nitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])[#6]")

            reactant_has_unsaturated = any(
                mol and mol.HasSubstructMatch(unsaturated_nitrile_pattern) for mol in reactants
            )
            product_has_saturated = product and product.HasSubstructMatch(saturated_nitrile_pattern)

            if reactant_has_unsaturated and product_has_saturated:
                has_reduction_after_knoevenagel = True
                print(f"Found reduction of α,β-unsaturated nitrile at depth {depth}")
                continue

        # Check for alkylation at geminal dinitrile carbon after reduction
        if has_reduction_after_knoevenagel and not has_alkylation_after_reduction and i > 1:
            geminal_dinitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])[#6]")
            alkylated_dinitrile_pattern = Chem.MolFromSmarts("[#6]([C]#[N])([C]#[N])([#6])[#6]")

            reactant_has_geminal = any(
                mol and mol.HasSubstructMatch(geminal_dinitrile_pattern) for mol in reactants
            )
            product_has_alkylated = product and product.HasSubstructMatch(
                alkylated_dinitrile_pattern
            )

            if reactant_has_geminal and product_has_alkylated:
                has_alkylation_after_reduction = True
                print(f"Found alkylation at geminal dinitrile carbon at depth {depth}")
                continue

    # Check if the complete strategy was found
    strategy_detected = (
        has_knoevenagel and has_reduction_after_knoevenagel and has_alkylation_after_reduction
    )

    if strategy_detected:
        print("Detected complete malononitrile-based linear synthesis strategy")
    else:
        print("Did not detect complete malononitrile-based linear synthesis strategy")
        print(
            f"Knoevenagel: {has_knoevenagel}, Reduction: {has_reduction_after_knoevenagel}, Alkylation: {has_alkylation_after_reduction}"
        )

    return strategy_detected
