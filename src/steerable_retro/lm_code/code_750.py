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
    This function detects if the synthesis involves a sequence of aromatic functionalizations.
    """
    # Track consecutive aromatic modifications
    consecutive_aromatic_mods = 0
    max_consecutive_mods = 0

    def dfs_traverse(node):
        nonlocal consecutive_aromatic_mods, max_consecutive_mods

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an aromatic functionalization
                # Simple heuristic: both reactants and products contain aromatic rings
                # and there's a change in functional groups
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and any(reactant_mols):
                    aromatic_pattern = Chem.MolFromSmarts("c")

                    reactants_aromatic = any(
                        r and r.HasSubstructMatch(aromatic_pattern) for r in reactant_mols
                    )
                    product_aromatic = product_mol.HasSubstructMatch(aromatic_pattern)

                    if reactants_aromatic and product_aromatic:
                        # This is a simplification - in a real implementation, you'd want to check
                        # that the same aromatic ring is being modified
                        consecutive_aromatic_mods += 1
                    else:
                        max_consecutive_mods = max(max_consecutive_mods, consecutive_aromatic_mods)
                        consecutive_aromatic_mods = 0

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    max_consecutive_mods = max(max_consecutive_mods, consecutive_aromatic_mods)

    # Consider it a sequence if there are 3 or more consecutive modifications
    has_aromatic_sequence = max_consecutive_mods >= 3
    if has_aromatic_sequence:
        print(
            f"Aromatic functionalization sequence detected with {max_consecutive_mods} consecutive modifications"
        )

    return has_aromatic_sequence
