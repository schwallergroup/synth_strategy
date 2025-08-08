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
    Detects if the synthesis involves thiazole ring formation from a thiourea and a bromomethyl ketone.
    """
    thiazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal thiazole_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain thiourea and bromomethyl ketone patterns
                thiourea_pattern = Chem.MolFromSmarts("[NH2][C](=[S])[NH]")
                bromomethyl_ketone_pattern = Chem.MolFromSmarts("[C](=[O])[CH2][Br]")

                has_thiourea = False
                has_bromomethyl_ketone = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(thiourea_pattern):
                            has_thiourea = True
                        if mol.HasSubstructMatch(bromomethyl_ketone_pattern):
                            has_bromomethyl_ketone = True

                # Check if product contains thiazole
                product_mol = Chem.MolFromSmiles(product)
                thiazole_pattern = Chem.MolFromSmarts("c1scnc1")

                if (
                    product_mol
                    and has_thiourea
                    and has_bromomethyl_ketone
                    and product_mol.HasSubstructMatch(thiazole_pattern)
                ):
                    print("Detected thiazole ring formation from thiourea and bromomethyl ketone")
                    thiazole_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiazole_formation_detected
