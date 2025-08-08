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
    Detects if the synthetic route uses a late-stage amide coupling strategy
    where a carboxylic acid and amine are joined in one of the final steps.
    """
    final_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal final_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid pattern in reactants
                acid_pattern = Chem.MolFromSmarts("[C,c](=O)[OH]")
                # Check for amine pattern in reactants
                amine_pattern = Chem.MolFromSmarts("[NH2][c,C]")
                # Check for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[C,c](=O)[NH][c,C]")

                acid_present = False
                amine_present = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(acid_pattern):
                            acid_present = True
                        if mol.HasSubstructMatch(amine_pattern):
                            amine_present = True

                product_mol = Chem.MolFromSmiles(product)
                amide_present = False
                if product_mol:
                    amide_present = product_mol.HasSubstructMatch(amide_pattern)

                if acid_present and amine_present and amide_present:
                    print(f"Found late-stage amide coupling at depth {depth}")
                    final_amide_coupling = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return final_amide_coupling
