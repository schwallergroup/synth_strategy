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
    Detects if the synthesis route uses a protection-deprotection strategy for carboxylic acids.
    Specifically, looks for carboxylic acid → ester → carboxylic acid sequence.
    """
    # Track protection and deprotection events
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for esterification (protection)
            acid_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8;H,-]")
            ester_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8][#6]")

            # Check for protection: acid → ester
            reactant_has_acid = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(acid_pattern):
                    reactant_has_acid = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            if reactant_has_acid and product_mol and product_mol.HasSubstructMatch(ester_pattern):
                protection_events.append(depth)
                print(f"Detected carboxylic acid protection at depth {depth}")

            # Check for deprotection: ester → acid
            reactant_has_ester = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(ester_pattern):
                    reactant_has_ester = True
                    break

            if reactant_has_ester and product_mol and product_mol.HasSubstructMatch(acid_pattern):
                deprotection_events.append(depth)
                print(f"Detected carboxylic acid deprotection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    if protection_events and deprotection_events:
        # Ensure protection happens before deprotection in the synthesis direction
        # (which means higher depth number for protection in retrosynthesis)
        for prot_depth in protection_events:
            for deprot_depth in deprotection_events:
                if prot_depth > deprot_depth:
                    print("Confirmed protection-deprotection sequence for carboxylic acid")
                    return True

    return False
