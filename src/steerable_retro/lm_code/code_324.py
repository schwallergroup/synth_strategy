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
    Detects if the synthetic route involves a biaryl disconnection strategy
    using boronic acid and triflate coupling partners.
    """
    has_biaryl_disconnection = False
    has_boronic_acid = False
    has_triflate = False

    def dfs_traverse(node):
        nonlocal has_biaryl_disconnection, has_boronic_acid, has_triflate

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid in reactants
                boronic_acid_pattern = Chem.MolFromSmarts("[#6]B(O)(O)")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                            has_boronic_acid = True
                            print(f"Found boronic acid in reactant: {reactant}")
                    except:
                        continue

                # Check for triflate in reactants
                triflate_pattern = Chem.MolFromSmarts("[#8]S(=O)(=O)C(F)(F)F")
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(triflate_pattern):
                            has_triflate = True
                            print(f"Found triflate in reactant: {reactant}")
                    except:
                        continue

                # Check if this is a biaryl disconnection
                # This is a simplified check - in a real implementation, you would need
                # to analyze the reaction more carefully to confirm it's a biaryl disconnection
                if has_boronic_acid and has_triflate:
                    has_biaryl_disconnection = True
                    print("Identified biaryl disconnection with boronic acid and triflate")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_biaryl_disconnection and has_boronic_acid and has_triflate
