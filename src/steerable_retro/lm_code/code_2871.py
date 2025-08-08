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
    This function detects a strategy involving sequential protection and deprotection
    of carboxylic acids.
    """
    # Track protection and deprotection events by depth
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if not product_mol:
                    return

                # Check for protection (carboxylic acid → ester)
                carboxylic_pattern = Chem.MolFromSmarts("C(=O)[OH]")
                tbutyl_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)")
                methyl_pattern = Chem.MolFromSmarts("COC(=O)")

                # Protection: reactant has COOH, product has ester
                if any(m and m.HasSubstructMatch(carboxylic_pattern) for m in reactant_mols if m):
                    if product_mol.HasSubstructMatch(
                        tbutyl_pattern
                    ) or product_mol.HasSubstructMatch(methyl_pattern):
                        protection_events.append(depth)
                        print(f"Found protection event at depth {depth}")

                # Deprotection: reactant has ester, product has COOH
                if any(
                    m
                    and (m.HasSubstructMatch(tbutyl_pattern) or m.HasSubstructMatch(methyl_pattern))
                    for m in reactant_mols
                    if m
                ):
                    if product_mol.HasSubstructMatch(carboxylic_pattern):
                        deprotection_events.append(depth)
                        print(f"Found deprotection event at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    has_protection = len(protection_events) > 0
    has_deprotection = len(deprotection_events) > 0

    # Check if protection happens before deprotection (higher depth = earlier in synthesis)
    sequential_strategy = False
    if has_protection and has_deprotection:
        # Check if at least one protection event happens at a higher depth than a deprotection event
        if max(protection_events) > min(deprotection_events):
            sequential_strategy = True

    print(f"Protection events at depths: {protection_events}")
    print(f"Deprotection events at depths: {deprotection_events}")
    print(f"Sequential protection-deprotection strategy detected: {sequential_strategy}")

    return sequential_strategy
