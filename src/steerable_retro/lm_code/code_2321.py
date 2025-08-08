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
    Detects a synthetic strategy involving biaryl formation via cross-coupling.
    """
    has_biaryl_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_biaryl_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for biaryl formation (two aromatic rings connected)
            # Look for reactions where separate aromatic rings in reactants become connected in product
            if len(reactants) >= 2:
                # Check if product contains biaryl motif
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        biaryl_pattern = Chem.MolFromSmarts(
                            "[c]!@[c]"
                        )  # Non-fused aromatic C connected to another aromatic C
                        if prod_mol.HasSubstructMatch(biaryl_pattern):
                            # Check if reactants have aromatic rings
                            aromatic_reactants = 0
                            for r in reactants:
                                r_mol = Chem.MolFromSmiles(r)
                                if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("c")):
                                    aromatic_reactants += 1

                            if aromatic_reactants >= 2:
                                has_biaryl_formation = True
                                print(f"Found biaryl formation at depth {depth}")
                except:
                    pass  # Handle parsing errors gracefully

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_biaryl_formation
