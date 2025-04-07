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
    Detects if the synthesis route involves a protection/deprotection sequence for an alcohol.
    Specifically looks for TBDMS (tert-butyldimethylsilyl) protection.
    """
    protection_events = []

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and all(r for r in reactant_mols):
                # TBDMS pattern
                tbdms_pattern = Chem.MolFromSmarts("[#14]([#6])([#6])[#8]")

                # Check for protection (TBDMS addition)
                if product_mol.HasSubstructMatch(tbdms_pattern) and not any(
                    r.HasSubstructMatch(tbdms_pattern) for r in reactant_mols
                ):
                    protection_events.append(("protection", depth))
                    print(f"TBDMS protection detected at depth {depth}")

                # Check for deprotection (TBDMS removal)
                if any(
                    r.HasSubstructMatch(tbdms_pattern) for r in reactant_mols
                ) and not product_mol.HasSubstructMatch(tbdms_pattern):
                    protection_events.append(("deprotection", depth))
                    print(f"TBDMS deprotection detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    has_protection = any(event == "protection" for event, _ in protection_events)
    has_deprotection = any(event == "deprotection" for event, _ in protection_events)

    if has_protection and has_deprotection:
        print("Complete alcohol protection/deprotection sequence detected")

    return has_protection and has_deprotection
