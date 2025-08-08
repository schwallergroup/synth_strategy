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
    This function detects a sequential functional group interconversion pathway:
    nitrile → carboxylic acid → ester
    """
    # Track reactions in sequence
    nitrile_to_acid_found = False
    acid_to_ester_found = False

    def dfs_traverse(node):
        nonlocal nitrile_to_acid_found, acid_to_ester_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile to acid conversion
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                    acid_pattern = Chem.MolFromSmarts("[C](=[O])[O;H1]")
                    ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

                    # Check for nitrile to acid conversion
                    if any(
                        r and r.HasSubstructMatch(nitrile_pattern) for r in reactant_mols
                    ) and product_mol.HasSubstructMatch(acid_pattern):
                        nitrile_to_acid_found = True
                        print("Found nitrile to carboxylic acid conversion")

                    # Check for acid to ester conversion
                    if any(
                        r and r.HasSubstructMatch(acid_pattern) for r in reactant_mols
                    ) and product_mol.HasSubstructMatch(ester_pattern):
                        acid_to_ester_found = True
                        print("Found carboxylic acid to ester conversion")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitrile_to_acid_found and acid_to_ester_found
