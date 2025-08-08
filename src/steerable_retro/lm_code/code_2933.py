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
    Detects if the synthetic route involves nucleophilic aromatic substitution
    with thiocyanate as a leaving group.
    """
    found_snar_with_thiocyanate = False

    def dfs_traverse(node):
        nonlocal found_snar_with_thiocyanate

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains thiocyanate group
            thiocyanate_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]#[#7]")
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(thiocyanate_pattern):
                        # Check if the product doesn't have the thiocyanate or has it in a different position
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            # If product doesn't have thiocyanate or has fewer thiocyanate groups
                            if not prod_mol.HasSubstructMatch(thiocyanate_pattern) or len(
                                prod_mol.GetSubstructMatches(thiocyanate_pattern)
                            ) < len(mol.GetSubstructMatches(thiocyanate_pattern)):
                                # Check if the reaction involves an aromatic ring
                                aromatic_pattern = Chem.MolFromSmarts("a")
                                if mol.HasSubstructMatch(aromatic_pattern):
                                    print("Found SNAr with thiocyanate leaving group")
                                    found_snar_with_thiocyanate = True
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_snar_with_thiocyanate
