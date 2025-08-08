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
    This function detects a synthetic strategy involving nitro-olefin intermediates
    with C-C bond formations adjacent to the nitro group.
    """
    nitro_olefin_present = False
    cc_bond_formation_with_nitro = False

    def dfs_traverse(node):
        nonlocal nitro_olefin_present, cc_bond_formation_with_nitro

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            products = [Chem.MolFromSmiles(p) for p in products_part.split(".") if p]

            # Check for nitro-olefin pattern in reactants or products
            nitro_olefin_pattern = Chem.MolFromSmarts("C=C[N+](=[O])[O-]")
            for mol in reactants + products:
                if mol and mol.HasSubstructMatch(nitro_olefin_pattern):
                    print("Found nitro-olefin pattern")
                    nitro_olefin_present = True

            # Check for C-C bond formation adjacent to nitro group
            if len(products) == 1 and len(reactants) >= 1:
                for r in reactants:
                    if r and r.HasSubstructMatch(Chem.MolFromSmarts("[N+](=[O])[O-]")):
                        print("Found potential C-C bond formation with nitro group")
                        cc_bond_formation_with_nitro = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = nitro_olefin_present and cc_bond_formation_with_nitro
    print(f"Nitro-olefin strategy detected: {result}")
    return result
