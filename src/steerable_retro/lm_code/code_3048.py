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
    Detects if the synthesis uses a strategy of sequential C-C bond formations,
    specifically looking for biaryl coupling followed by olefin formation.
    """
    # Track C-C bond formations and their sequence
    cc_bond_formations = []

    def dfs_traverse(node, depth=0):
        nonlocal cc_bond_formations

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product is not None:
                # Check for biaryl coupling
                biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                if product.HasSubstructMatch(biaryl_pattern):
                    # Check if this bond is being formed
                    if not all(
                        r is not None and r.HasSubstructMatch(biaryl_pattern) for r in reactants
                    ):
                        cc_bond_formations.append(("biaryl", depth))
                        print(f"Found biaryl C-C bond formation at depth {depth}")

                # Check for olefin formation
                alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")
                if product.HasSubstructMatch(alkene_pattern):
                    # Check if this bond is being formed
                    if not all(
                        r is not None and r.HasSubstructMatch(alkene_pattern) for r in reactants
                    ):
                        cc_bond_formations.append(("olefin", depth))
                        print(f"Found olefin C=C bond formation at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check for the specific sequence: biaryl coupling followed by olefin formation
    has_biaryl = any(bond_type == "biaryl" for bond_type, _ in cc_bond_formations)
    has_olefin = any(bond_type == "olefin" for bond_type, _ in cc_bond_formations)

    # Check if biaryl occurs at a higher depth (earlier in synthesis) than olefin
    correct_sequence = False
    if has_biaryl and has_olefin:
        biaryl_depth = max(
            depth for bond_type, depth in cc_bond_formations if bond_type == "biaryl"
        )
        olefin_depth = max(
            depth for bond_type, depth in cc_bond_formations if bond_type == "olefin"
        )
        correct_sequence = biaryl_depth > olefin_depth

    strategy_present = has_biaryl and has_olefin and correct_sequence

    if strategy_present:
        print("Detected sequential C-C bond formation strategy (biaryl followed by olefin)")
    else:
        print("Sequential C-C bond formation strategy not detected")

    return strategy_present
