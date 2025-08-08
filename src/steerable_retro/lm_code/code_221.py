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
    This function detects if the synthesis involves formation of a biaryl ether bond,
    likely via Suzuki-Miyaura coupling.
    """
    has_biaryl_ether_formation = False

    def dfs_traverse(node):
        nonlocal has_biaryl_ether_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant has a boronic acid and another has a hydroxyl group
            boronic_acid_pattern = Chem.MolFromSmarts("[c]B(O)O")
            hydroxyl_pattern = Chem.MolFromSmarts("[c][O][H]")

            # Check if product has biaryl ether linkage
            biaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")

            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(biaryl_ether_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(
                            boronic_acid_pattern
                        ) or reactant_mol.HasSubstructMatch(hydroxyl_pattern):
                            has_biaryl_ether_formation = True
                            print("Detected biaryl ether formation via Suzuki coupling")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_biaryl_ether_formation
