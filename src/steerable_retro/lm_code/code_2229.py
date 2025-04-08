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
    This function detects a linear synthesis route featuring fluorinated aromatics,
    Weinreb amide as a ketone precursor, and late-stage Grignard addition.
    """
    has_fluorinated_aromatic = False
    has_weinreb_amide = False
    has_grignard_reaction = False
    has_late_stage_grignard = False

    def dfs_traverse(node, depth=0):
        nonlocal has_fluorinated_aromatic, has_weinreb_amide, has_grignard_reaction, has_late_stage_grignard

        if node["type"] == "mol":
            # Check for fluorinated aromatic in molecules
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for fluorinated aromatic
                    fluorinated_aromatic_pattern = Chem.MolFromSmarts("[c][F]")
                    if mol.HasSubstructMatch(fluorinated_aromatic_pattern):
                        has_fluorinated_aromatic = True
                        print(f"Found fluorinated aromatic in molecule: {node['smiles']}")

                    # Check for Weinreb amide
                    weinreb_amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]([#6])[#8][#6]")
                    if mol.HasSubstructMatch(weinreb_amide_pattern):
                        has_weinreb_amide = True
                        print(f"Found Weinreb amide in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Grignard reaction (contains Mg in reactants)
                if any("Mg" in reactant for reactant in reactants):
                    has_grignard_reaction = True
                    print(f"Found Grignard reaction: {rsmi}")

                    # Check if it's a late-stage Grignard (depth <= 1)
                    if depth <= 1:
                        has_late_stage_grignard = True
                        print(f"Found late-stage Grignard reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if all key features are present
    strategy_present = (
        has_fluorinated_aromatic
        and has_weinreb_amide
        and has_grignard_reaction
        and has_late_stage_grignard
    )
    print(f"Strategy detection result: {strategy_present}")
    return strategy_present
