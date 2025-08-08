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
    This function detects thiazole ring formation via Hantzsch thiazole synthesis
    from an α-bromoketone and thiourea/thioamide.
    """
    thiazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal thiazole_formation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains thiazole
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    thiazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#6][#16]1")
                    if product_mol.HasSubstructMatch(thiazole_pattern):

                        # Check if reactants contain α-bromoketone and thiourea/thioamide
                        bromo_ketone_pattern = Chem.MolFromSmarts("[#6][#6](=[#8])[#6][Br]")
                        thiourea_pattern = Chem.MolFromSmarts("[#7][#6]([#7])=[#16]")

                        has_bromo_ketone = False
                        has_thiourea = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                if reactant_mol.HasSubstructMatch(bromo_ketone_pattern):
                                    has_bromo_ketone = True
                                if reactant_mol.HasSubstructMatch(thiourea_pattern):
                                    has_thiourea = True

                        if has_bromo_ketone and has_thiourea:
                            print("Detected Hantzsch thiazole synthesis")
                            thiazole_formation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return thiazole_formation_detected
