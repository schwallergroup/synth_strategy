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
    Detects a specific functional group transformation sequence:
    COOH → COOMe → CH2OH → CH2OMs → CH2N
    """
    # Track if we've found each transformation in the sequence
    found_acid_to_ester = False
    found_ester_to_alcohol = False
    found_alcohol_to_mesylate = False
    found_mesylate_to_amine = False

    def dfs_traverse(node):
        nonlocal found_acid_to_ester, found_ester_to_alcohol, found_alcohol_to_mesylate, found_mesylate_to_amine

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and all(r for r in reactants):
                    # Check for carboxylic acid to ester transformation
                    acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8;H1]")
                    ester_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8]-[#6]")

                    if any(
                        r.HasSubstructMatch(acid_pattern) for r in reactants
                    ) and product.HasSubstructMatch(ester_pattern):
                        print("Found carboxylic acid to ester transformation")
                        found_acid_to_ester = True

                    # Check for ester to alcohol transformation
                    if any(r.HasSubstructMatch(ester_pattern) for r in reactants):
                        alcohol_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#8;H1]")
                        if product.HasSubstructMatch(alcohol_pattern):
                            print("Found ester to alcohol transformation")
                            found_ester_to_alcohol = True

                    # Check for alcohol to mesylate transformation
                    alcohol_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#8;H1]")
                    mesylate_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#8]-[#16](=[#8])(=[#8])-[#6]")

                    if any(
                        r.HasSubstructMatch(alcohol_pattern) for r in reactants
                    ) and product.HasSubstructMatch(mesylate_pattern):
                        print("Found alcohol to mesylate transformation")
                        found_alcohol_to_mesylate = True

                    # Check for mesylate to amine transformation
                    if any(r.HasSubstructMatch(mesylate_pattern) for r in reactants):
                        amine_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#7](-[#6])-[#6]")
                        if product.HasSubstructMatch(amine_pattern):
                            print("Found mesylate to amine transformation")
                            found_mesylate_to_amine = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # The sequence is present if we found all the transformations
    sequence_present = (
        found_acid_to_ester
        and found_ester_to_alcohol
        and found_alcohol_to_mesylate
        and found_mesylate_to_amine
    )
    print(f"Functional group transformation sequence detected: {sequence_present}")
    return sequence_present
