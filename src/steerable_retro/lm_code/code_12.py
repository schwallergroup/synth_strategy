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
    This function detects a synthetic strategy involving phenol alkylation
    to form aryl ethers.
    """
    phenol_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal phenol_alkylation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for phenol and aryl ether
                phenol_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]")
                aryl_ether_pattern = Chem.MolFromSmarts(
                    "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]-[#6]"
                )
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")

                # Check reactants and products
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                # Check if phenol and alkyl halide are reactants and aryl ether is product
                if product_mol and product_mol.HasSubstructMatch(aryl_ether_pattern):
                    if any(
                        mol and mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols
                    ) and any(
                        mol and mol.HasSubstructMatch(alkyl_halide_pattern) for mol in reactant_mols
                    ):
                        print("Detected phenol alkylation")
                        phenol_alkylation_detected = True

                # Also check reverse (retrosynthetic direction)
                if product_mol and any(
                    mol and mol.HasSubstructMatch(aryl_ether_pattern) for mol in reactant_mols
                ):
                    if product_mol.HasSubstructMatch(phenol_pattern):
                        print("Detected aryl ether cleavage (reverse of phenol alkylation)")
                        phenol_alkylation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return phenol_alkylation_detected
