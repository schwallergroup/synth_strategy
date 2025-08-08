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
    Detects a synthetic strategy involving introduction of a CF3-containing aryl group.
    """
    has_cf3_aryl_introduction = False

    def dfs_traverse(node):
        nonlocal has_cf3_aryl_introduction

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for CF3-aryl group in product
            cf3_aryl_pattern = Chem.MolFromSmarts("[c]-[C]([F])([F])[F]")
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and product_mol.HasSubstructMatch(cf3_aryl_pattern):
                # Check if CF3-aryl group is in reactants
                cf3_in_reactants = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(cf3_aryl_pattern):
                        cf3_in_reactants = True
                        break

                if cf3_in_reactants:
                    has_cf3_aryl_introduction = True
                    print("Detected CF3-aryl group introduction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_cf3_aryl_introduction
