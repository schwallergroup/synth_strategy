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


def main(route, min_interconversions=2):
    """
    Detects if the synthesis route includes multiple interconversions between
    carboxylic acids and esters (protection/deprotection)
    """
    interconversion_count = 0

    def dfs_traverse(node):
        nonlocal interconversion_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit mol objects
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and all(reactant_mols):
                    # SMARTS for carboxylic acid
                    acid_pattern = Chem.MolFromSmarts("[#6][C](=[O])[OH]")

                    # SMARTS for ester
                    ester_pattern = Chem.MolFromSmarts("[#6][C](=[O])[O][#6]")

                    # Check for acid to ester conversion
                    has_acid_in_reactants = any(
                        r.HasSubstructMatch(acid_pattern) for r in reactant_mols
                    )
                    has_ester_in_product = product_mol.HasSubstructMatch(ester_pattern)

                    # Check for ester to acid conversion
                    has_ester_in_reactants = any(
                        r.HasSubstructMatch(ester_pattern) for r in reactant_mols
                    )
                    has_acid_in_product = product_mol.HasSubstructMatch(acid_pattern)

                    if (has_acid_in_reactants and has_ester_in_product) or (
                        has_ester_in_reactants and has_acid_in_product
                    ):
                        interconversion_count += 1
                        print(f"Found acid/ester interconversion, total: {interconversion_count}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return interconversion_count >= min_interconversions
