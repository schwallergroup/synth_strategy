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
    This function detects if the route uses a nitrile -> amine -> amide
    functional group interconversion sequence.
    """
    has_nitrile_to_amine = False
    has_amine_to_amide = False

    def dfs_traverse(node):
        nonlocal has_nitrile_to_amine, has_amine_to_amide

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            amine_pattern = Chem.MolFromSmarts("[N;H2][C]")
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

            product_mol = Chem.MolFromSmiles(product)

            # Check for nitrile -> amine conversion
            if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                nitrile_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitrile_pattern):
                        nitrile_in_reactants = True
                        break

                if nitrile_in_reactants:
                    has_nitrile_to_amine = True
                    print("Nitrile to amine conversion detected")

            # Check for amine -> amide conversion
            if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                amine_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(amine_pattern):
                        amine_in_reactants = True
                        break

                if amine_in_reactants:
                    has_amine_to_amide = True
                    print("Amine to amide conversion detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    result = has_nitrile_to_amine and has_amine_to_amide
    print(f"Nitrile -> amine -> amide sequence: {result}")
    return result
