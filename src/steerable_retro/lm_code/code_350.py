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
    Detects if the synthesis route uses a halogenation followed by coupling strategy.
    Specifically, looks for bromination followed by a C-C bond formation.
    """
    # Track halogenation and coupling events with their depths
    halogenation_events = []
    coupling_events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for halogenation (specifically bromination)
            if "Br" in product and any("Br" in reactant for reactant in reactants):
                # Simple check for bromination - product has bromine and at least one reactant has bromine
                # This suggests bromine transfer to the substrate
                halogenation_events.append(depth)
                print(f"Detected potential bromination at depth {depth}")

            # Check for coupling reactions (C-C bond formation)
            # For Suzuki coupling, look for boron-containing reactant
            boron_pattern = Chem.MolFromSmarts("[#5]")
            has_boron = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(boron_pattern):
                    has_boron = True
                    break

            # If we have a boron-containing reactant, it's likely a coupling reaction
            if has_boron:
                coupling_events.append(depth)
                print(f"Detected potential coupling reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both halogenation and coupling events
    if halogenation_events and coupling_events:
        # Ensure halogenation happens before coupling in the synthesis direction
        # (which means higher depth number for halogenation in retrosynthesis)
        for hal_depth in halogenation_events:
            for coup_depth in coupling_events:
                if hal_depth > coup_depth:
                    print("Confirmed halogenation-coupling sequence")
                    return True

    return False
