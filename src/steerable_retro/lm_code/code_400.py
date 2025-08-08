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
    This function detects if the synthesis route uses a Heck coupling reaction
    for C-C bond formation between aryl and alkene components.
    """
    has_heck_coupling = False

    def dfs_traverse(node):
        nonlocal has_heck_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Heck coupling pattern
            reactants_mol = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and len(reactants_mol) >= 2:
                # Check if one reactant has aryl group
                has_aryl = any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts("c")) for mol in reactants_mol if mol
                )

                # Check if one reactant has alkene
                has_alkene = any(
                    mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C]"))
                    for mol in reactants_mol
                    if mol
                )

                # Check if product has aryl-alkene connection
                has_aryl_alkene = product_mol.HasSubstructMatch(Chem.MolFromSmarts("c/[C]=[C]"))

                if has_aryl and has_alkene and has_aryl_alkene:
                    has_heck_coupling = True
                    print("Detected Heck coupling for C-C bond formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_heck_coupling
