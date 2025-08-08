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
    Detects if the synthetic route uses early benzylation to connect aromatic fragments
    via C-C bond formation
    """
    found_benzylation = False
    early_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal found_benzylation, early_stage

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Look for patterns indicating benzylation
                aryl_pattern = Chem.MolFromSmarts("c")
                benzyl_pattern = Chem.MolFromSmarts("c[CX4]")

                has_aryl = False
                has_benzyl = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(aryl_pattern):
                            has_aryl = True
                        if mol.HasSubstructMatch(benzyl_pattern):
                            has_benzyl = True

                # Check if this might be a benzylation reaction
                if has_aryl and (
                    has_benzyl or "Zn" in rsmi
                ):  # Zinc often indicates benzylation conditions
                    found_benzylation = True

                    # Check if this is in the first half of synthesis (depth > 3)
                    if depth > 3:
                        early_stage = True
                        print(f"Found early-stage benzylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_benzylation and early_stage
