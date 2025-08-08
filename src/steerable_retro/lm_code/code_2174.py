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
    Detects a synthetic strategy involving an oxidation step (thioether to sulfone)
    followed by a reduction step (nitro to amine)
    """
    # Initialize tracking variables
    has_oxidation = False
    oxidation_depth = -1
    has_reduction = False
    reduction_depth = -1

    # SMARTS patterns
    thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-c")
    sulfone_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])(=[#8])-c")
    nitro_pattern = Chem.MolFromSmarts("[$([#7](=[#8])~[#8])]")
    amine_pattern = Chem.MolFromSmarts("[NH2]-c")

    def dfs_traverse(node, depth=0):
        nonlocal has_oxidation, oxidation_depth, has_reduction, reduction_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants = parts[0].split(".")
                    product = parts[2]

                    # Convert to RDKit molecules
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    # Check for oxidation (thioether to sulfone)
                    if product_mol and product_mol.HasSubstructMatch(sulfone_pattern):
                        if any(
                            r and r.HasSubstructMatch(thioether_pattern) for r in reactant_mols if r
                        ):
                            has_oxidation = True
                            oxidation_depth = depth
                            print(f"Detected oxidation (thioether to sulfone) at depth {depth}")

                    # Check for reduction (nitro to amine)
                    if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                        if any(
                            r and r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r
                        ):
                            has_reduction = True
                            reduction_depth = depth
                            print(f"Detected reduction (nitro to amine) at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present (oxidation followed by reduction)
    strategy_present = (
        has_oxidation and has_reduction and oxidation_depth > reduction_depth
    )  # Higher depth means earlier in synthesis

    if strategy_present:
        print("Detected oxidation-reduction sequence strategy")

    return strategy_present
