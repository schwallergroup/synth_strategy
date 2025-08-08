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
    Detects a linear synthesis strategy where each step primarily involves
    functional group transformations rather than carbon skeleton changes.
    """
    carbon_skeleton_changes = 0
    functional_group_changes = 0

    def dfs_traverse(node):
        nonlocal carbon_skeleton_changes, functional_group_changes

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is primarily a functional group transformation
                # by comparing carbon counts
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product)

                if all(reactant_mols) and product_mol:
                    # Count carbons in main reactant (assuming first reactant is main)
                    main_reactant_carbons = sum(
                        1 for atom in reactant_mols[0].GetAtoms() if atom.GetAtomicNum() == 6
                    )
                    product_carbons = sum(
                        1 for atom in product_mol.GetAtoms() if atom.GetAtomicNum() == 6
                    )

                    # If carbon count is similar, it's likely a functional group transformation
                    if abs(main_reactant_carbons - product_carbons) <= 1:
                        functional_group_changes += 1
                    else:
                        carbon_skeleton_changes += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is primarily functional group transformations
    total_reactions = carbon_skeleton_changes + functional_group_changes
    is_linear_fg_strategy = (
        total_reactions > 0 and functional_group_changes / total_reactions >= 0.66
    )

    print(f"Linear functional group transformation strategy: {is_linear_fg_strategy}")
    print(f"Functional group transformations: {functional_group_changes}")
    print(f"Carbon skeleton changes: {carbon_skeleton_changes}")

    return is_linear_fg_strategy
