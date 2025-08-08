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
    This function detects N-alkylation of indole with a benzyl halide.
    """
    n_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for indole NH and benzyl halide patterns
                indole_nh_pattern = Chem.MolFromSmarts("[nH]1c2ccccc2cc1")
                benzyl_halide_pattern = Chem.MolFromSmarts("c[CH2][Br,Cl,I,F]")
                n_alkylated_indole_pattern = Chem.MolFromSmarts("[n]([CH2]c)1c2ccccc2cc1")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(m and m.HasSubstructMatch(indole_nh_pattern) for m in reactant_mols)
                    and any(m and m.HasSubstructMatch(benzyl_halide_pattern) for m in reactant_mols)
                ):
                    if product_mol.HasSubstructMatch(n_alkylated_indole_pattern):
                        print("Indole N-alkylation detected")
                        n_alkylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_alkylation_detected
