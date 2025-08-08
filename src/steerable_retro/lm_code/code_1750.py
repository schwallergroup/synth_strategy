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
    This function detects nitro group reduction to amine.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

            for reactant in reactants:
                try:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol and r_mol.HasSubstructMatch(nitro_pattern):
                        # Check if product has amine where nitro was
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol:
                            # Count nitro groups in reactant
                            nitro_count_r = len(r_mol.GetSubstructMatches(nitro_pattern))
                            # Count nitro groups in product
                            nitro_count_p = len(p_mol.GetSubstructMatches(nitro_pattern))

                            # Check if nitro groups decreased and amines increased
                            amine_pattern = Chem.MolFromSmarts("[NH2]")
                            amine_count_r = len(r_mol.GetSubstructMatches(amine_pattern))
                            amine_count_p = len(p_mol.GetSubstructMatches(amine_pattern))

                            if nitro_count_p < nitro_count_r and amine_count_p > amine_count_r:
                                print("Found nitro reduction to amine")
                                nitro_reduction_found = True
                except:
                    continue

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_found
