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
    Detects a synthetic strategy involving multiple silyl protections and
    phosphate ester formation.
    """
    # Track if we found the key transformations
    found_silyl_protection = False
    found_phosphate_formation = False

    def dfs_traverse(node):
        nonlocal found_silyl_protection, found_phosphate_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not product or not all(reactants):
                return

            # Check for silyl protection
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH]")
            silyl_ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[Si]")
            tbdms_pattern = Chem.MolFromSmarts("[#6]-[#8]-[Si]([#6])([#6])[C]([#6])([#6])[#6]")

            for reactant in reactants:
                if reactant.HasSubstructMatch(alcohol_pattern):
                    if product and (
                        product.HasSubstructMatch(silyl_ether_pattern)
                        or product.HasSubstructMatch(tbdms_pattern)
                    ):
                        print("Found silyl protection")
                        found_silyl_protection = True

            # Check for phosphate ester formation
            phosphate_pattern = Chem.MolFromSmarts("[#6]-[#8]-P(=O)(-[#8]-[#6])-[#8]-[#6]")
            if product and product.HasSubstructMatch(phosphate_pattern):
                print("Found phosphate ester formation")
                found_phosphate_formation = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both key transformations were found
    return found_silyl_protection and found_phosphate_formation
