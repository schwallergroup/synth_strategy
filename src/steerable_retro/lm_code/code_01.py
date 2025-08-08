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
    This function detects a synthetic strategy involving multiple SNAr reactions.
    SNAr reactions typically involve displacement of a halogen on an electron-deficient
    aromatic ring by a nitrogen nucleophile.
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr pattern: halogen on aromatic being replaced by nitrogen
                # Look for nitrogen nucleophile in reactants
                nitrogen_nucleophile = False
                for reactant in reactants:
                    if "[NH" in reactant or "[N:" in reactant:
                        nitrogen_nucleophile = True
                        break

                # Look for halogen on aromatic in reactants
                halogen_aromatic = False
                for reactant in reactants:
                    if ("[Cl" in reactant or "[Br" in reactant or "[F" in reactant) and (
                        "[c" in reactant or "[n" in reactant
                    ):
                        halogen_aromatic = True
                        break

                # Check if C-N bond formed where halogen was
                if nitrogen_nucleophile and halogen_aromatic:
                    # This is a simplified check - a more robust implementation would use actual reaction mapping
                    snar_count += 1
                    print(f"SNAr reaction detected at depth: {node.get('depth', 'unknown')}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple SNAr reactions are detected
    return snar_count >= 2
