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
    Detects a synthetic strategy involving multiple protection-deprotection cycles
    for amine functional groups using different protecting groups.
    """
    # Track protection-deprotection cycles
    protection_deprotection_cycles = []

    # Define patterns for common amine protecting groups
    protecting_groups = {
        "trifluoroacetamide": Chem.MolFromSmarts("[N]-C(=O)-C([F])([F])([F])"),
        "carbamate": Chem.MolFromSmarts("[N]-C(=O)-O-[C]"),
        "amide": Chem.MolFromSmarts("[N]-C(=O)-[C,c]"),
        "boc": Chem.MolFromSmarts("[N]-C(=O)-O-C(C)(C)C"),
        "cbz": Chem.MolFromSmarts("[N]-C(=O)-O-C-c"),
    }

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            product_mol = Chem.MolFromSmiles(product_part)

            # Check for protection (adding protecting group)
            for pg_name, pg_pattern in protecting_groups.items():
                if product_mol and product_mol.HasSubstructMatch(pg_pattern):
                    # Check if any reactant doesn't have the protecting group
                    has_reactant_without_pg = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and not reactant_mol.HasSubstructMatch(pg_pattern):
                            has_reactant_without_pg = True
                            break

                    if has_reactant_without_pg:
                        protection_deprotection_cycles.append(
                            {"type": "protection", "group": pg_name, "depth": depth}
                        )
                        print(f"Found {pg_name} protection at depth {depth}")

            # Check for deprotection (removing protecting group)
            for pg_name, pg_pattern in protecting_groups.items():
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(pg_pattern):
                        if product_mol and not product_mol.HasSubstructMatch(pg_pattern):
                            protection_deprotection_cycles.append(
                                {"type": "deprotection", "group": pg_name, "depth": depth}
                            )
                            print(f"Found {pg_name} deprotection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze protection-deprotection cycles
    complete_cycles = {}
    for event in protection_deprotection_cycles:
        pg = event["group"]
        if pg not in complete_cycles:
            complete_cycles[pg] = {"protection": False, "deprotection": False}

        if event["type"] == "protection":
            complete_cycles[pg]["protection"] = True
        elif event["type"] == "deprotection":
            complete_cycles[pg]["deprotection"] = True

    # Count complete cycles (both protection and deprotection)
    cycle_count = sum(
        1 for pg, cycle in complete_cycles.items() if cycle["protection"] and cycle["deprotection"]
    )

    # Strategy requires at least 2 complete protection-deprotection cycles
    result = cycle_count >= 2
    print(f"Multiple amine protection-deprotection cycles: {result} (found {cycle_count} cycles)")
    return result
