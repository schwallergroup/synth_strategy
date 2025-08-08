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
    Detects N-methylation reaction in the synthetic route.
    """
    found_n_methylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_methylation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Pattern for secondary amine (NH)
            sec_amine_pattern = Chem.MolFromSmarts("[#6]-[NH]-[#6]")
            # Pattern for tertiary amine with methyl (N-CH3)
            tert_amine_pattern = Chem.MolFromSmarts("[#6]-[N]([CH3])-[#6]")

            has_sec_amine = False

            # Check reactants for secondary amine
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(sec_amine_pattern):
                        has_sec_amine = True
                        break
                except:
                    continue

            # Check product for N-methylated amine
            if has_sec_amine:
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(tert_amine_pattern):
                        found_n_methylation = True
                        print(f"Found N-methylation at depth {depth}")
                except:
                    pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_n_methylation
