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
    This function detects a synthetic strategy involving chain extension via
    C-C bond formation using an activated methylene (e.g., chloromethyl).
    """
    chain_extension_detected = False

    def dfs_traverse(node):
        nonlocal chain_extension_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for chain extension via alkylation
            # Look for activated methylene (chloromethyl) in reactants
            chloromethyl_pattern = Chem.MolFromSmarts("[Cl][CH2][#6]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and any(
                r and r.HasSubstructMatch(chloromethyl_pattern) for r in reactant_mols
            ):
                # Check if the product has a new C-C bond that wasn't in the reactants
                # This is a simplification - in a real implementation, you'd need more sophisticated
                # atom mapping and comparison
                if len(reactants) >= 2:  # Need at least two reactants for C-C bond formation
                    print(f"Detected potential chain extension via alkylation: {rsmi}")
                    chain_extension_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Chain extension via alkylation strategy detected: {chain_extension_detected}")
    return chain_extension_detected
