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
    This function detects a strategy involving C-C bond formations adjacent to
    electron-withdrawing groups (nitro, carbonyl).
    """
    cc_bond_formations = 0

    def dfs_traverse(node):
        nonlocal cc_bond_formations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Check for C-C bond formation adjacent to EWG
            # EWG patterns: nitro, carbonyl
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
            carbonyl_pattern = Chem.MolFromSmarts("C(=O)")

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            products = [Chem.MolFromSmiles(p) for p in products_part.split(".") if p]

            # This is a simplified check - in a real implementation, you would need to
            # analyze the actual reaction mechanism to detect bond formations
            for p in products:
                if p:
                    if p.HasSubstructMatch(nitro_pattern) and any(
                        r.HasSubstructMatch(nitro_pattern) for r in reactants if r
                    ):
                        print("Potential C-C bond formation near nitro group")
                        cc_bond_formations += 1

                    if p.HasSubstructMatch(carbonyl_pattern) and any(
                        r.HasSubstructMatch(carbonyl_pattern) for r in reactants if r
                    ):
                        print("Potential C-C bond formation near carbonyl group")
                        cc_bond_formations += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = cc_bond_formations >= 2
    print(f"C-C bond formation adjacent to EWG strategy detected: {result}")
    return result
