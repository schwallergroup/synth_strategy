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
    Detects if the synthesis includes a Boc protection step.
    Looks for formation of a tert-butoxycarbonyl group on nitrogen.
    """
    boc_protection_found = False

    def dfs_traverse(node):
        nonlocal boc_protection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains Boc group
            boc_pattern = Chem.MolFromSmarts("[#7]-[C](=[O])-[O]-C([C])([C])[C]")

            # Check if reactants include Boc anhydride or similar
            boc_reagent_pattern = re.compile(r"CC\(C\)\(C\)OC\(=O\)O")

            boc_reagent_found = False
            for reactant in reactants:
                if boc_reagent_pattern.search(reactant):
                    boc_reagent_found = True
                    break

            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(boc_pattern) and boc_reagent_found:
                    boc_protection_found = True
                    print("Detected Boc protection step")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_protection_found
